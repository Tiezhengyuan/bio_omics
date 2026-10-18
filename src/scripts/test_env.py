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
    import pyrosetta
    pyrosetta.init("-mute all")
    print('\n\nPyRosetta version:', pyrosetta.version())


if __name__ == "__main__":
    main()