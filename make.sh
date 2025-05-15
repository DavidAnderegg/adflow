export PETSC_ARCH=$PETSC_ARCH_COMPLEX
make -f Makefile_CS


export PETSC_ARCH=$PETSC_ARCH_REAL
make


