try:
    import cantera as ct
    from sdtoolbox.cp import cpsolve # noqa
except ImportError:
    print("Something wrong with the import")

print("Using Cantera version: ", ct.__version__)
