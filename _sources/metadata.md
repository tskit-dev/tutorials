---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{currentmodule} tskit
```


(sec_tutorial_metadata)=

# Working with Metadata

Metadata is information associated with entities that {program}`tskit` doesn't use or
interpret, but which is useful to pass on to downstream analysis such as sample ids,
dates etc. (see {ref}`sec_metadata` for a full discussion). Each
{ref}`table<sec_tables_api_table>` has a {class}`MetadataSchema` which details the
contents and encoding of the metadata for each row. A metadata schema is a JSON document
that conforms to [JSON Schema](https://json-schema.org/understanding-json-schema/)
(The full schema for tskit is at {ref}`sec_metadata_schema_schema`). Here we use an
{ref}`example tree sequence<sec_intro_downloading_datafiles>`
which contains some demonstration metadata:

```{code-cell} ipython3
:tags: [remove-cell]
import msprime
import tskit

def metadata():
    tables = msprime.sim_ancestry(4).dump_tables()
    tables.individuals.metadata_schema = tskit.MetadataSchema(
    {'additionalProperties': False,
     'codec': 'json',
     'properties': {'accession': {'description': 'ENA accession number',
                                  'type': 'string'},
                    'pcr': {'description': 'Was PCR used on this sample',
                            'name': 'PCR Used',
                            'type': 'boolean'}},
     'required': ['accession', 'pcr'],
     'type': 'object'}
    )
    md = [
        {'accession': 'ERS0001', 'pcr': True},
        {'accession': 'ERS0002', 'pcr': True},
        {'accession': 'ERS0003', 'pcr': True},
        {'accession': 'ERS0004', 'pcr': False},
    ]
    table = tables.individuals
    copy = table.copy()
    table.clear()
    for m, row in zip(md, copy):
        table.append(row.replace(metadata=m))
    ts = tables.tree_sequence()
    ts.dump("data/metadata.trees")

def create_notebook_data():
    metadata()

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```


```{code-cell} ipython3
import tskit
import json

ts = tskit.load("data/metadata.trees")
```
(sec_tutorial_metadata_reading)=

## Reading metadata and schemas

Metadata is automatically decoded using the schema when accessed via a
{class}`TreeSequence` or {class}`TableCollection` Python API. For example:

```{code-cell} ipython3
print("Metadata for individual 0:", ts.individual(0).metadata)  # Tree sequence access
print("Metadata for individual 0:", ts.tables.individuals[0].metadata)  # Table access
```

Viewing the {class}`MetadataSchema` for a table can help with understanding
its metadata, as it can contain descriptions and constraints:

```{code-cell} ipython3
ts.table_metadata_schemas.individual
```

The same schema can be accessed via a {attr}`~IndividualTable.metadata_schema` attribute
on each table (printed prettily here using ``json.dumps``)
 
```{code-cell} ipython3
schema = ts.tables.individuals.metadata_schema
print(json.dumps(schema.asdict(), indent=4))  # Print with indentations
```

The top-level metadata and schemas for the entire tree sequence are similarly
accessed with {attr}`TreeSequence.metadata` and {attr}`TreeSequence.metadata_schema`.

:::{note}
If there is no schema (i.e. it is equal to ``MetadataSchema(None)``) for a table
or top-level metadata, then no decoding is performed and ``bytes`` will be returned.
:::

(sec_tutorial_metadata_modifying)=

## Modifying metadata and schemas

If you are creating or modifying a tree sequence by changing the underlying tables,
you may want to record or add to the metadata. If the change fits into the same schema,
this is relatively simple, you can follow the
{ref}`description of minor table edits<sec_tables_editing_minor>` in the
{ref}`sec_tables` tutorial. However if it requires a change to the schema, this must be
done first, as it is then used to validate and encode the metadata. 

Schemas in tskit are held in a {class}`MetadataSchema`.
A Python dict representation of the schema is passed to its constructor, which
will validate the schema. Here are a few examples: the first one allows arbitrary fields
to be added, the second one (which will construct the schema we printed above) does not:

```{code-cell} ipython3
basic_schema = tskit.MetadataSchema({'codec': 'json'})

complex_schema = tskit.MetadataSchema({
    'codec': 'json',
    'additionalProperties': False,
    'properties': {'accession': {'description': 'ENA accession number',
                                 'type': 'string'},
                   'pcr': {'description': 'Was PCR used on this sample',
                           'name': 'PCR Used',
                           'type': 'boolean'}},
    'required': ['accession', 'pcr'],
    'type': 'object',
})
```

This {class}`MetadataSchema` can then be assigned to a table or the top-level
tree sequence e.g. {attr}`~IndividualTable.metadata_schema`:

```{code-cell} ipython3
tables = tskit.TableCollection(sequence_length=1)  # make a new, empty set of tables
tables.individuals.metadata_schema = complex_schema
```

This will overwrite any existing schema. Note that this will not validate any existing
metadata against the new schema. Now that the table has a schema, calls to
{meth}`~IndividualTable.add_row` will validate and encode the metadata:

```{code-cell} ipython3
row_id = tables.individuals.add_row(0, metadata={"accession": "Bob1234", "pcr": True})
print(f"Row {row_id} added to the individuals table")
```

If we try to add metadata that doesn't fit the schema, such as accidentally using a
string instead of a proper Python boolean, we'll get an error:

```{code-cell} ipython3
:tags: [raises-exception, output_scroll]
tables.individuals.add_row(0, metadata={"accession": "Bob1234", "pcr": "false"})
```

and because we set ``additionalProperties`` to ``False`` in the schema, an error is
also raised if we attempt to add new fields:

```{code-cell} ipython3
:tags: [raises-exception, output_scroll]
tables.individuals.add_row(0, metadata={"accession": "Bob1234", "pcr": True, "newKey": 25})
```


To set the top-level metadata, just assign it. Validation and encoding happen as
specified by the top-level metadata schema

```{code-cell} ipython3
tables.metadata_schema = basic_schema  # Allows new fields to be added that are not validated
tables.metadata = {"mean_coverage": 200.5}
print(tables.metadata)
```

:::{note}
*Provenance* information, detailing the origin of the data, modification timestamps,
and (ideally) how the tree sequence can be reconstructed, should go in
{ref}`sec_provenance`, not metadata.
:::

To modify a schema --- for example to add a key --- first get the dict representation,
modify, then write back:

```{code-cell} ipython3
schema_dict = tables.individuals.metadata_schema.schema
schema_dict["properties"]["newKey"] = {"type": "integer"}
tables.individuals.metadata_schema = tskit.MetadataSchema(schema_dict)
# Now this will work:
new_id = tables.individuals.add_row(metadata={'accession': 'abc123', 'pcr': False, 'newKey': 25})
print(tables.individuals[new_id].metadata)
```

To modify the metadata of rows in tables use the {ref}`sec_tutorial_metadata_bulk`.

(sec_tutorial_metadata_viewing_raw)=

## Viewing raw metadata

If you need to see the raw (i.e. bytes) metadata, you just need to remove the
schema, for instance:

```{code-cell} ipython3
individual_table = tables.individuals.copy()  # don't change the original tables.individual

print("Metadata:\n", individual_table[0].metadata)

individual_table.metadata_schema = tskit.MetadataSchema(None)
print("\nRaw metadata:\n", individual_table[0].metadata)
```

(sec_tutorial_metadata_bulk)=

## Metadata for bulk table methods

In the interests of efficiency each table's {meth}`~NodeTable.packset_metadata` method,
as well as the more general {meth}`~NodeTable.set_columns` and
{meth}`~NodeTable.append_columns` methods, do not attempt to validate or encode metadata.
You can call {meth}`MetadataSchema.validate_and_encode_row` directly to prepare metadata
for these methods:

```{code-cell} ipython3
metadata_column = [
    {"accession": "etho1234", "pcr": True},
    {"accession": "richard1235", "pcr": False},
    {"accession": "albert1236", "pcr": True},
]
encoded_metadata_column = [
    tables.individuals.metadata_schema.validate_and_encode_row(r) for r in metadata_column
]
md, md_offset = tskit.pack_bytes(encoded_metadata_column)
tables.individuals.set_columns(flags=[0, 0, 0], metadata=md, metadata_offset=md_offset)
tables.individuals
```

Or if all columns do not need to be set:

```{code-cell} ipython3
tables.individuals.packset_metadata(
    [tables.individuals.metadata_schema.validate_and_encode_row(r) for r in metadata_column]
)
```

(sec_tutorial_metadata_binary)=

## Binary metadata

To disable the validation and encoding of metadata and store raw bytes pass ``None`` to
{class}`MetadataSchema`

```{code-cell} ipython3
tables.populations.metadata_schema = tskit.MetadataSchema(None)
tables.populations.add_row(metadata=b"SOME CUSTOM BYTES #!@")
print(tables.populations[0].metadata)
```
