# These commands will be executed after doing "make install"
message("Running PostInstall.cmake...")

# Create a pkg-config file to allow other programs to link with MAST more easily
set(DEST_DIR "${CMAKE_INSTALL_PREFIX}")
set(DEST_LIB_DIR ${CMAKE_INSTALL_LIBDIR})
set(PRIVATE_LIBS "${PRIVATE_LIBS} -lmast${CMAKE_POSTFIX}") # For pkg-config file generation

CONFIGURE_FILE(
    "${ROOT_DIR}/mast.pc.in" 
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig/mast${CMAKE_POSTFIX}.pc"
    @ONLY)

# Create a symbolic link to the libmast_postfix.so named libmast.so for legacy support
execute_process(COMMAND ln -s "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libmast${CMAKE_POSTFIX}.so" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libmast.so")

# Recommend some post-install steps
message("\nIt is recommended to take the following steps:")
message("\tAdd \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\" to LD_LIBRARY_PATH")
message("\tAdd \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig\" to PKG_CONFIG_PATH")
message("\tAdd \"${CMAKE_INSTALL_PREFIX}/include\" to CPATH")

message("\nNOTE: Using pkg-config, the flags to link against MAST can be obtained by running \"pkg-config --cflags --libs mast${CMAKE_POSTFIX}\"\n")

message("Completed running PostInstall.cmake")
