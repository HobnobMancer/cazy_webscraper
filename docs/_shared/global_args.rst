.. note::
   ``--version``, ``--citation``, ``-l``/``--log``, ``--sql_echo`` and ``-v``/``--verbose`` belong
   to ``cazy_webscraper`` itself, not to any subcommand, so they must be given **before** the
   subcommand name:

   .. code-block:: bash

      cazy_webscraper -v -l run.log get_uniprot_data cazy/cazyme.db

   Placing them after the subcommand is an error. See :doc:`migration` if you are converting a
   version 2 command.
