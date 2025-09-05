# preset that turns on just a few packages for UEFEX


set(ALL_PACKAGES MOLECULE UEF UEFEX)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(LAMMPS_SIZES "bigbig" CACHE STRING "" FORCE)
