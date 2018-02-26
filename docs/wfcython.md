
# A Wright-Fisher simulation implemented in C via Cython.

This tutorial implements a Wright-Fisher simulation with mutation and recombination using [Cython](http://www.cython.org).  Cython is two things:

* A grammer/dialect of Python that allows static typic of C/C++ types.
* A static compiler to turn the Cython grammer in to C or C++ code to compile into a Python extension module.

Cython has a learning curve of its own. A lot of what is shown below reflects best practices.  For those, we refer you to the Cython documentation.

Here, we avoid all use of [numpy](http://www.numpy.org) until we have to talk to [msprime](http://msprime.readthedocs.io).  We replace all numpy functionality with the equivalent routines from the excellent GNU Scientific Library, or [GSL](https://www.gnu.org/software/gsl/doc/html/index.html). Yes, numpy is fast!  Numpy is written in C!  But, numpy has to talk back and forth to Python, meaning we can out-perform it by writing routines that execute completely on the C side.

This example is closer to reality for those working in lower-level languages.  First, we must build our world, which means defining data types (structs, in this case), functions acting on those types, and a bunch of auxillary code to manage memory and handle errors.  After all that, we can code up the `simplify` and `evolve` functions. Such is the price of speed.

First, we load an extension allowing us to write Cython in a notebook:


```python
%load_ext Cython
```

The next code block is our Cython code.  The notebook environment will magically compile it from Cython to C, from C to a compiled Python module, and load the module into memory.

The code block is unavoidably long.  It defines a single Python function called `evolve`, which may then be run later in the notebook.  All of the functions marked `cdef` are visible only as C functions by other functions within the module.  This limited function scope is why we must write everything in one giant block.

Some comments:

* We use [CythonGSL](https://github.com/twiecki/CythonGSL) to get access to GSL types and functions in Cython.  Conda or pip install it if you want to use it for your projects.
* The memory managment and error handling on the C side is minimal, and the error handling in particular is naive.
* The recombination function is 100% executed in C.  Mutation is very close to that, except for where we use the `dict` to manage the infinitely-many sites mutation model.
* We get copy-free transfer (I think...) from C to numpy arrays via Cython's typed memory views.
* The `cimport` commands below bring names into scope. It is considered best practice to only `cimport` the symbols you use, but that quickly gets tedious here, and I gave myself a break and imported everything from `gsl_vector`.


```cython
%%cython -3 -lgsl -lgslcblas -lm

import msprime
import numpy as np
from collections import namedtuple
import pickle
cimport numpy as np
from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, realloc, free
from libc.stdint cimport int32_t, uint32_t

from cython_gsl.gsl_rng cimport gsl_rng
from cython_gsl.gsl_rng cimport gsl_rng_mt19937
from cython_gsl.gsl_rng cimport gsl_rng_alloc
from cython_gsl.gsl_rng cimport gsl_rng_free
from cython_gsl.gsl_rng cimport gsl_rng_set
from cython_gsl.gsl_rng cimport gsl_rng_uniform
from cython_gsl.gsl_random cimport gsl_ran_flat
from cython_gsl.gsl_random cimport gsl_ran_poisson
from cython_gsl.gsl_vector cimport *
from cython_gsl.gsl_sort cimport gsl_sort_vector

cdef int32_t * malloc_int32_t(size_t n):
    return <int32_t*>malloc(n*sizeof(int32_t))

cdef int32_t * realloc_int32_t(void * x, size_t n):
    return <int32_t*>realloc(x,n*sizeof(int32_t))

cdef double * malloc_double(size_t n):
    return <double*>malloc(n*sizeof(double))

cdef double * realloc_double(double * x, size_t n):
    return <double*>realloc(<double *>x,n*sizeof(double))

cdef struct Mutations:
    double * pos
    int32_t * time
    int32_t * node
    size_t next_mutation, capacity
    
cdef int init_Mutations(Mutations * m):
    m.next_mutation = 0
    m.capacity = 10000
    m.pos = malloc_double(m.capacity)
    if m.pos == NULL:
        return -1
    m.time = malloc_int32_t(m.capacity)
    if m.time == NULL:
        return -1
    m.node = malloc_int32_t(m.capacity)
    if m.node == NULL:
        return -1
    return 0

cdef int realloc_Mutations(Mutations * m):
    m.capacity *= 2
    m.pos = realloc_double(m.pos,
                          m.capacity)
    if m.pos == NULL:
        return -1
    m.time = realloc_int32_t(m.time,
                            m.capacity)
    if m.time == NULL:
        return -1
    m.node = realloc_int32_t(m.node,
                            m.capacity)
    if m.node == NULL:
        return -1
    return 0

cdef void free_Mutations(Mutations * m):
    free(m.pos)
    free(m.time)
    free(m.node)
    m.next_mutation = 0
    m.capacity = 10000
    
cdef int add_mutation(double pos,
                     int32_t generation,
                     int32_t node,
                     Mutations * m):
    cdef int rv = 0
    if m.next_mutation+1 >= m.capacity:
        rv = realloc_Mutations(m)
        if rv != 0:
            return rv
    m.pos[m.next_mutation] = pos
    m.time[m.next_mutation] = generation
    m.node[m.next_mutation] = node
    m.next_mutation+=1
    return rv
    
cdef struct Nodes:
    double * time
    size_t next_node, capacity
    
cdef int init_Nodes(Nodes * n):
    n.next_node = 0
    n.capacity = 10000
    n.time = malloc_double(n.capacity)
    if n.time == NULL:
        return -1
    return 0

cdef int realloc_Nodes(Nodes * n):
    n.capacity *= 2
    n.time = realloc_double(n.time,
                            n.capacity)
    if n.time == NULL:
        return -1
    return 0
    
cdef void free_Nodes(Nodes * n):
    if n.time != NULL:
        free(n.time)
    n.next_node = 0
    n.capacity = 10000

cdef int add_node(double t, Nodes *n):
    cdef int rv = 0
    if n.next_node >= n.capacity:
        rv = realloc_Nodes(n)
        if rv != 0:
            return rv
    n.time[n.next_node] = t
    n.next_node+=1
    return rv
    
cdef struct Edges:
    double *left
    double *right
    int32_t *parent
    int32_t *child
    size_t next_edge, capacity
    
cdef int init_Edges(Edges * e):
    e.next_edge = 0
    e.capacity = 10000
    e.left = malloc_double(e.capacity)
    if e.left == NULL:
        return -1
    e.right = malloc_double(e.capacity)
    if e.right == NULL:
        return -1
    e.parent = malloc_int32_t(e.capacity)
    if e.parent == NULL:
        return -1
    e.child = malloc_int32_t(e.capacity)
    if e.child == NULL:
        return -1
    return 0
   
cdef int realloc_Edges(Edges * e):
    e.capacity *= 2
    e.left = realloc_double(e.left,e.capacity)
    if e.left == NULL:
        return -1
    e.right = realloc_double(e.right,e.capacity)
    if e.right == NULL:
        return -1
    e.parent = realloc_int32_t(e.parent,e.capacity)
    if e.parent == NULL:
        return -1
    e.child = realloc_int32_t(e.child,e.capacity)
    if e.child == NULL:
        return -1
    return 0

cdef void free_Edges(Edges * e):
    free(e.left)
    free(e.right)
    free(e.parent)
    free(e.child)
    e.next_edge = 0
    e.capacity = 10000
    
cdef int add_edge(double left, double right,
             int32_t parent, int32_t child,
             Edges * edges):
    cdef int rv=0
    if edges.next_edge+1 >= edges.capacity:
        rv = realloc_Edges(edges)
        if rv != 0:
            return rv
        
    edges.left[edges.next_edge] = left
    edges.right[edges.next_edge] = right
    edges.parent[edges.next_edge] = parent
    edges.child[edges.next_edge] = child
    edges.next_edge += 1
    return rv

cdef struct Tables:
    Nodes nodes
    Edges edges
    Mutations mutations
    gsl_rng * rng
    
cdef int init_Tables(Tables * t, int seed):
    cdef int rv = 0
    rv = init_Nodes(&t.nodes)
    if rv != 0:
        return rv
    rv = init_Edges(&t.edges)
    if rv != 0:
        return rv
    rv = init_Mutations(&t.mutations)
    if rv != 0:
        return rv
    t.rng = gsl_rng_alloc(gsl_rng_mt19937)
    if t.rng == NULL:
        return -1
    gsl_rng_set(t.rng, seed)
    return rv

cdef void free_Tables(Tables * t):
    free_Nodes(&t.nodes)
    free_Edges(&t.edges)
    free_Mutations(&t.mutations)
    gsl_rng_free(t.rng)
    
cdef int infsites(double mu, int32_t generation,
                  int32_t next_offspring_index,
                  Tables * tables,
                  dict lookup):
    cdef unsigned nmut = gsl_ran_poisson(tables.rng, mu)
    cdef unsigned i = 0
    cdef double pos
    cdef int rv = 0
    for i in range(nmut):
        pos = gsl_rng_uniform(tables.rng)
        while pos in lookup:
            pos = gsl_rng_uniform(tables.rng)
        rv = add_mutation(pos,
                         generation,
                         next_offspring_index,
                         &tables.mutations)
        if rv != 0:
            return rv
        lookup[pos] = True
    return rv

cdef int value_present_vector(gsl_vector * v, double x,
                              size_t start, size_t stop):
    cdef size_t i
    for i in range(start,stop):
        if gsl_vector_get(v,i)==x:
            return 1
    return 0

cdef int poisson_recombination(double r,
                               size_t pg1, size_t pg2,
                               int32_t next_offspring_id,
                               Tables * tables):
    cdef unsigned nbreaks = gsl_ran_poisson(tables.rng, r)
    cdef gsl_vector * b = NULL
    cdef unsigned i = 0,drew_zero=0
    cdef double x
    cdef int rv = 0
    cdef double left,right
    cdef int32_t p
    if nbreaks == 0:
        # The parent passes the entire region onto the child
        rv = add_edge(0.0,1.0,pg1,
                      next_offspring_id,&tables.edges)
        if rv != 0:
            return rv
    else:
        b = gsl_vector_calloc(nbreaks+2) 
        while i < nbreaks:
            x = gsl_rng_uniform(tables.rng)
            while value_present_vector(b,x,1,i+1)==1:
                x = gsl_rng_uniform(tables.rng)
            gsl_vector_set(b,i+1,x)
            if x == 0:
                drew_zero=1
            i += 1
        gsl_vector_set(b,b.size-1,1.0)
        gsl_sort_vector(b)
        if drew_zero == 1:
            pg1,pg2 = pg2,pg1
        for i in range(b.size-1):
            left = gsl_vector_get(b,i)
            right = gsl_vector_get(b,i+1)
            rv = add_edge(left,right,pg1,
                          next_offspring_id,&tables.edges)
            if rv != 0:
                gsl_vector_free(b)
                return rv
            pg1,pg2 = pg2,pg1
    gsl_vector_free(b)
    return 0

cdef int make_offspring(double mu, double r,
                        size_t generation,
                        size_t pg1, size_t pg2,
                        int32_t next_offspring_index,
                        dict lookup,
                        Tables * tables):
    cdef int rv
    rv = poisson_recombination(r,pg1,pg2,
                               next_offspring_index,
                               tables)
    if rv != 0:
        return -2
                
    rv = infsites(mu,generation+1,
                  next_offspring_index,
                  tables,lookup)
    if rv != 0:
        return -3
            
    rv = add_node(generation+1, &tables.nodes)
    if rv != 0:
        return -4
   
    return 0

cdef void handle_error_code(int error, Tables * tables):
    """
    Only to be called after make_offspring
    """
    if error == 0:
        return
    print("Error occurred")
    free_Tables(tables)
    if error == -2:
        raise RuntimeError("error during recombination")
    elif error == -2:
        raise RuntimeError("error during mutation")
    elif error == -4:
        raise RuntimeError("erorr adding nodes")
    else:
        raise ValueError("invalid error code")

MutationMetadata = namedtuple('MutationMetadata',['pos','origin'])

cdef int simplify(Tables * tables, 
            double dt,
            object nodes,
            object edges,
            object sites,
            object mutations):
    cdef int rv = 0,gap
    
    cdef size_t i
    cdef np.ndarray[double,ndim=1] dview,lview,rview
    cdef np.ndarray[int32_t,ndim=1] pview,cview
    # Reverse time for our new nodes
    cdef gsl_vector_view vt
    vt = gsl_vector_view_array(tables.nodes.time,<size_t>tables.nodes.next_node)
    cdef double tmax,tmin
    gsl_vector_minmax(&vt.vector,&tmin,&tmax)
    for i in range(tables.nodes.next_node):
        tables.nodes.time[i] = -1.0*(tables.nodes.time[i]-tmax)
    gsl_vector_minmax(&vt.vector,&tmin,&tmax)
    
    nodes.set_columns(time=nodes.time+dt,flags=nodes.flags)
    dview = np.asarray(<double[:tables.nodes.next_node]>tables.nodes.time)
    gap=nodes.time.min()-tmax
    if gap != 1:
        return -1
    nodes.append_columns(time=dview,
                        flags=np.ones(tables.nodes.next_node,dtype=np.uint32))
    
    
    lview = np.asarray(<double[:tables.edges.next_edge]>tables.edges.left)
    rview = np.asarray(<double[:tables.edges.next_edge]>tables.edges.right)
    pview = np.asarray(<int32_t[:tables.edges.next_edge]>tables.edges.parent)
    cview = np.asarray(<int32_t[:tables.edges.next_edge]>tables.edges.child)
    edges.append_columns(left=lview,
                        right=rview,
                        parent=pview,
                        child=cview)
    
    # We are trying to be as fast as possible here,
    # so we'll use the more cumbersome 
    # append_columns interface instead of the 
    # much slower (but easier to understand)
    # add_rows
    cdef size_t nsites = len(sites)
    dview = np.asarray(<double[:tables.mutations.next_mutation]>tables.mutations.pos)
    sites.append_columns(position=dview,
                        ancestral_state=np.zeros(len(dview),dtype=np.int8)+ord('0'),
                        ancestral_state_offset=np.arange(len(dview)+1,dtype=np.uint32))
    pview = np.asarray(<int32_t[:tables.mutations.next_mutation]>tables.mutations.node)
    mutations.append_columns(site=np.arange(tables.mutations.next_mutation,
                                            dtype=np.int32)+nsites,
                            node=pview,
                            derived_state=np.ones(len(dview),
                                                  dtype=np.int8)+ord('0'),
                            derived_state_offset=np.arange(len(dview)+1,
                                                          dtype=np.uint32)
                            )
    
    msprime.sort_tables(nodes=nodes,edges=edges,
                       sites=sites,mutations=mutations)
    
    samples = np.where(nodes.time == 0)[0]
    msprime.simplify_tables(samples=samples.tolist(),
                            nodes=nodes,
                           edges=edges,
                            sites=sites,
                            mutations=mutations)
    #print(len(sites),len(mutations))
    
    # "clear" our temp containers
    tables.nodes.next_node = 0
    tables.mutations.next_mutation = 0
    tables.edges.next_edge = 0
                          
    return rv

def evolve(int N, int ngens, double theta, double rho, int gc, int seed):
    nodes = msprime.NodeTable()
    edges = msprime.EdgeTable()
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    
    cdef double mu = theta/<double>(4*N)
    cdef double r = rho/<double>(4*N)
    
    cdef int rv
    cdef size_t i, generation
    cdef Tables tables
    rv = init_Tables(&tables, seed)
    if rv != 0:
        free_Tables(&tables)
        raise RuntimeError("could not initialize tables")
        
    for i in range(2*<size_t>N):
        nodes.add_row(time=0.0,
                      flags=msprime.NODE_IS_SAMPLE)
        
    cdef int32_t next_offspring_index, first_parental_index
    next_offspring_index = len(nodes)
    first_parental_index = 0
    cdef size_t parent1, parent2,pindex
    cdef int32_t p1g1, p1g2, p2g1, p2g2
    cdef dict lookup = {}
    cdef size_t last_gen_gc = 0
    for generation in range(<size_t>(ngens)):
        if generation>0 and generation%gc == 0.0:
            rv = simplify(&tables,
                         generation-last_gen_gc,
                         nodes,edges,sites,mutations)
            if rv != 0:
                free_Tables(&tables)
                raise RuntimeError("simplification error")
            # lookup = {i:True for i in sites.position}
            last_gen_gc=generation
            next_offspring_index = len(nodes)
            first_parental_index = 0
        else:
            first_parental_index = next_offspring_index - 2*N
            
        for pindex in range(0,2*N,2):
            parent1=<size_t>gsl_ran_flat(tables.rng,0.0,<double>N)
            parent2=<size_t>gsl_ran_flat(tables.rng,0.0,<double>N)
            p1g1 = first_parental_index + 2*parent1
            p1g2 = p1g1 + 1
            p2g1 = first_parental_index + 2*parent2
            p2g2 = p2g1 + 1
            
            if gsl_rng_uniform(tables.rng) < 0.5:
                p1g1, p1g2 = p1g2, p1g1
            if gsl_rng_uniform(tables.rng) < 0.5:
                p2g1, p2g2 = p2g2, p2g1
                
            rv = make_offspring(mu,r,generation,
                                p1g1,p1g2,
                                next_offspring_index,
                                lookup,
                                &tables)
            handle_error_code(rv,&tables)
            next_offspring_index+=1
            rv = make_offspring(mu,r,generation,
                                p2g1,p2g2,
                                next_offspring_index,
                                lookup,
                                &tables)
            handle_error_code(rv,&tables)
            next_offspring_index+=1
        
    if tables.nodes.next_node > 0:
        rv=simplify(&tables,
                   generation+1-last_gen_gc,
                    nodes,edges,sites,mutations
                   )
        if rv == -1:
            free_Tables(&tables)
            raise RuntimeError("simplification error")
    
    free_Tables(&tables)
    return msprime.load_tables(nodes=nodes,edges=edges,
                               sites=sites,mutations=mutations)
```

    /Users/kevin/anaconda3/lib/python3.5/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
      from ._conv import register_converters as _register_converters



```python
%%time
ts = evolve(100, 1000, 100.0, 100.0, 1, 42)
```

    CPU times: user 1.76 s, sys: 214 ms, total: 1.98 s
    Wall time: 1.98 s


Make sure that output is invariant to how often we simplify:


```python
for gc in range(10,1000,29):
    ts2 = evolve(100, 1000, 100.0, 100.0, gc, 42)
    assert(ts2.tables.nodes == ts.tables.nodes)
    assert(ts2.tables.edges == ts.tables.edges)
    assert(ts2.tables.sites == ts.tables.sites)
    assert(ts2.tables.mutations == ts.tables.mutations)
```


```python
%%time
ts = evolve(1000,10000,1000.0,1000.0,1000,42)
```

    CPU times: user 16.8 s, sys: 2.15 s, total: 19 s
    Wall time: 19 s



```python
%%time
ts = evolve(10000,10000,1000.0,1000.0,1000,42)
```

    CPU times: user 1min 25s, sys: 12.9 s, total: 1min 38s
    Wall time: 1min 38s


# TODO

* Run neutral simulations, process with pylibseq, and compare output to msprime!


```python

```
