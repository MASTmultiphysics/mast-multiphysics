# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.example_driver.Debug:
PostBuild.mast.Debug: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/example_driver
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/example_driver:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/libmast.a\
	/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_dbg.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/example_driver


PostBuild.mast.Debug:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/libmast.a:
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/libmast.a


PostBuild.mast_tests.Debug:
PostBuild.mast.Debug: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/mast_tests
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/mast_tests:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/libmast.a\
	/opt/local/lib/libboost_unit_test_framework-mt.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_dbg.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/mast_tests


PostBuild.example_driver.Release:
PostBuild.mast.Release: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/example_driver
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/example_driver:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/libmast.a\
	/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/example_driver


PostBuild.mast.Release:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/libmast.a:
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/libmast.a


PostBuild.mast_tests.Release:
PostBuild.mast.Release: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/mast_tests
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/mast_tests:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/libmast.a\
	/opt/local/lib/libboost_unit_test_framework-mt.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/mast_tests


PostBuild.example_driver.MinSizeRel:
PostBuild.mast.MinSizeRel: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/example_driver
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/example_driver:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/libmast.a\
	/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/example_driver


PostBuild.mast.MinSizeRel:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/libmast.a:
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/libmast.a


PostBuild.mast_tests.MinSizeRel:
PostBuild.mast.MinSizeRel: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/mast_tests
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/mast_tests:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/libmast.a\
	/opt/local/lib/libboost_unit_test_framework-mt.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/mast_tests


PostBuild.example_driver.RelWithDebInfo:
PostBuild.mast.RelWithDebInfo: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/example_driver
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/example_driver:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/libmast.a\
	/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/example_driver


PostBuild.mast.RelWithDebInfo:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/libmast.a:
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/libmast.a


PostBuild.mast_tests.RelWithDebInfo:
PostBuild.mast.RelWithDebInfo: /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/mast_tests
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/mast_tests:\
	/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/libmast.a\
	/opt/local/lib/libboost_unit_test_framework-mt.dylib\
	/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib\
	/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib\
	/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib\
	/usr/lib/liblapack.dylib\
	/usr/lib/libblas.dylib\
	/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a\
	/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a\
	/opt/local/lib/libgcc/libgfortran.4.dylib\
	/opt/local/lib/openmpi-clang/libmpi_cxx.dylib\
	/opt/local/lib/openmpi-clang/libmpi.dylib
	/bin/rm -f /Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/mast_tests




# For each target create a dummy ruleso the target does not have to exist
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Debug/libmast.a:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/MinSizeRel/libmast.a:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/RelWithDebInfo/libmast.a:
/Users/manav/Documents/codes/inhouse/mast/mast_master/build/xcode/Release/libmast.a:
/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_dbg.dylib:
/Users/manav/Documents/codes/numerical/libmesh/lib/libmesh_opt.dylib:
/Users/manav/Documents/codes/numerical/optimization_codes/dot/build/libdot.a:
/Users/manav/Documents/codes/numerical/optimization_codes/gcmma/build/libgcmma.a:
/Users/manav/Documents/codes/numerical/petsc/lib/libpetsc.dylib:
/Users/manav/Documents/codes/numerical/slepc/lib/libslepc.dylib:
/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib:
/opt/local/lib/libboost_unit_test_framework-mt.dylib:
/opt/local/lib/libgcc/libgfortran.4.dylib:
/opt/local/lib/openmpi-clang/libmpi.dylib:
/opt/local/lib/openmpi-clang/libmpi_cxx.dylib:
/usr/lib/libblas.dylib:
/usr/lib/liblapack.dylib:
