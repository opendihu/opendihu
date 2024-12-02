Checkpointing
=============

OpenDiHu ships with a builtin Checkpointing framework with different
implementation backend that allow to switch OutputWriter and InputReader.

This framework is not yet enabled by default and is only enabled if the ``checkpointing`` object is present in the ``config``.

Configuration
_____________

.. code-block:: python

  config = {
      "checkpointing": {
            "directory": "state",,
            "interval": 10,
            "type": "hdf5-combined",
            "autoRestore": False,
            "checkpointToRestore": None,
       }
  }

.. list-table::
  :header-rows: 1

   * - variable
     - description
     - default
   * - directory
     - The directory where checkpoints will be stored
     - "state"
   * - interval
     - The timestep interval in which checkpoints will be created
     - 1
   * - type
     - The checkpointing backend which should be used. Available is "hdf5-combined", "hdf5-independent", "json-combined", "json-independent". The difference between combined and independent is that the combined checkpointing implementations combine results of all ranks into a single file, while the independent method writes each rank into an independent file.
     - "hdf5-combined"
   * - autoRestore
     - Checkpointing uses `SCR <https://scr.readthedocs.io/>`_ which has the feature to automatically restore the last checkpoint that was written. This option allows to toggle this feature.
     - False
   * - checkpointToRestore
     - Specify a specific checkpoint name to restore. The names are analog to the filenames written to the checkpointing directory.
     - None

Implementation
______________

The checkpointing implementation consistens of 3 Main components, the :doc:`output_writer`, input reader (currently only used for checkpointing) and in manager interface bringing both together.

For any output_writer to work it needs to be extended by a possibility to use the `uniqueName()` of the field variables, which reflects the full chain of nested solvers.
For reference see the HDF5 or JSON implementation.

Input reader also need to be implementated for each output writer which should be used.
For this are two different types for each output writer type (HDF5, JSON) of input readers in use, the partial and the full_dataset input readers.
The partial only reads the dataset for their own rank into a int or double vector, using the ``chunkDims`` attribute (see :doc:`output_writer`) so it should only be used in combination with the combined output writer.
The full_dataset input reader, will always read the full dataset into a vector.
