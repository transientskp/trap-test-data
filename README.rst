===================================
LOFAR Transients Pipeline Test Data
===================================

This repository contains the data required to run the test suite for the
`LOFAR Transients Pipeline (TraP) <https://github.com/transientskp/tkp/>`_.

Although available as a separate repository, it is available as a
``submodule`` with the TraP repository. To fetch a copy of the TraP and set up
the repository for testing::

    $ git clone git@github.com:transientskp/tkp.git
    $ cd tkp
    $ git submodule init
    $ git submodule update

Note that the submodule provides a pointer to a specific commit (SHA1) in the
test data repository. This pointer is versioned along with the TraP. In other
words, a given TraP version should always contain a pointer to the appropriate
set of test data for that particular version. If you change the TraP so that a
different set of data is needed, update the submodule appropriately. For
example::

    # Adding new test data
    $ cd tkp/tests/data
    # hack hack hack
    $ git commit -m "Updated test data"
    $ git push origin master

    # Update the reference in the TraP project
    $ cd ..
    $ git commit data -m "Updated test data"
    $ git push origin master
