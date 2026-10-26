# Python API

The package exports the optical model and the helpers below; everything else is
internal. Units are metres, nanoseconds and nanometres; frames are described in
the [user guide](USER_GUIDE.md#frames-and-units).

## Optical models

```{eval-rst}
.. autoclass:: obdeect.model.OpticalModel
   :members:

.. autofunction:: obdeect.model.trace_many

.. autoexception:: obdeect.formats.ModelError
```

## Analysis

```{eval-rst}
.. autofunction:: obdeect.psf.containment

.. autofunction:: obdeect.psf.summary
```

## Compilation

```{eval-rst}
.. autofunction:: obdeect.compiler.compile_model

.. autofunction:: obdeect.repository.read_telescope
```
