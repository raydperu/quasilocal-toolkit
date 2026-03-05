import qlt 
print("=== QuasiLocal Toolkit Python Module Test ===") 
print(f"Module version: {qlt.__version__}") 
print("\nMain module contents:") 
for item in dir(qlt): 
    if not item.startswith('_'): 
        print(f"  - {item}") 
 
print("\nExamining submodules:") 
submodules = ['axisymmetric_j', 'clebsch_komar', 'quasilocal_mass', 'surgical_flux', 'vorticity_bounds'] 
for sub in submodules: 
    if hasattr(qlt, sub): 
        submodule = getattr(qlt, sub) 
        print(f"\n{sub} contents:") 
        for attr in dir(submodule): 
            if not attr.startswith('_'): 
                print(f"    - {attr}") 
    else: 
        print(f"\n{sub} not found!") 
