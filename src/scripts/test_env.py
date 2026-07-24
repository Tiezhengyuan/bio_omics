#!/usr/bin/python
# python scripts/test_env.py

def main():
    #############
    import torch
    print("Torch:", torch.__version__)
    print("CUDA version:", torch.version.cuda)
    print("CUDA available:", torch.cuda.is_available())
    if torch.cuda.is_available():
        print("GPU:", torch.cuda.get_device_name(0))

    #############
    import pdbfixer
    print("\n\nPdbfixer version:", getattr(pdbfixer, "__version__", "unknown"))

    #############
    from openmm import version, Platform
    print("\n\nOpenMM version:", version.version)
    print("Available platforms:")
    for i in range(Platform.getNumPlatforms()):
        platform = Platform.getPlatform(i)
        print(f"  {i}: {platform.getName()}")
    try:
        cuda = Platform.getPlatformByName("CUDA")
        print("\n✓ CUDA platform is available")
    except Exception as e:
        print("\n✗ CUDA platform is NOT available")
        print(e)

    #############    
    import pyrosetta
    pyrosetta.init("-mute all")
    print('\n\nPyRosetta version:', pyrosetta.version())


if __name__ == "__main__":
    main()