
Source repository for the tskit tutorials site, 
[https://tskit-dev.github.io/tutorials](https://tskit-dev.github.io/tutorials/)

**UNDER CONSTRUCTION:** This is a very early version, and really just a way to 
explore some options for presenting this content.

## Organisation

The ``docs`` directory is a [GitHub pages](https://pages.github.com/) site. This 
means that all the Markdown files in this directory are automatically converted to 
HTML and made available on the website.
The source content for each 'chapter' is a Jupyter notebook in the ``notebooks``
directory. Notebooks are then converted to Markdown using ``nbconvert``, and placed 
in the ``docs`` directory.

## Converting a notebook

To convert a notebook to markdown, use the following:

```shell
$ jupyter nbconvert --to markdown --output-dir docs/ notebooks/NOTEBOOK_NAME.ipynb
```

When adding a new notebook to the site, you need to then add the files to 
git:

```shell
$ git add docs/NOTEBOOK_NAME*
```

Finally, update the ``docs/README.md`` to insert a link to the new page.

## TODO

- Need standardised titles including authorship. 
- Need some sort of citation mechanism. Perhaps [this](https://github.com/takluyver/cite2c)?
- Main page needs some content explaining what the site is for.
- Better template? We can use any Jekyll template, so it's quite flexible.
- It would also be nice to have a download link to the original notebook.

