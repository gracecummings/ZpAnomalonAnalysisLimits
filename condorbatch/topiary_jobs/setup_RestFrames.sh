# Source this script to set the environmental variables for this RestFrames installation
# 
# Users can automatically call this script at their shell login by adding to their
# .bashrc/.profile/.zlogin etc. file:
#    . $(restframes-config --prefix)/libexec/setup_RestFrames.sh
# assuming restframes-config is in their PATH, or
#    . PATH_OF_RESTFRAMES/libexec/setup_RestFrames.sh
#
# This script if for bash-like shells, see setup_RestFrames.csh for csh-like shells.

drop_from_path()
{
   # Assert that we got enough arguments
   if test $# -ne 2 ; then
      echo "drop_from_path: needs 2 arguments"
      return 1
   fi

   p=$1
   drop=$2

   newpath=`echo $p | sed -e "s;:${drop}:;:;g" \
                          -e "s;:${drop};;g"   \
                          -e "s;${drop}:;;g"   \
                          -e "s;${drop};;g"`
}

if [ -n "${RESTFRAMESSYS}" ] ; then
   old_restframessys=${RESTFRAMESSYS}
fi

RESTFRAMESSYS=$PWD/RestFrames

if [ -n "${old_restframessys}" ] ; then
   if [ -n "${PATH}" ]; then
      drop_from_path "$PATH" ${old_restframessys}/bin
      PATH=$newpath
   fi
   if [ -n "${LD_LIBRARY_PATH}" ]; then
      drop_from_path $LD_LIBRARY_PATH ${old_restframessys}/lib
      LD_LIBRARY_PATH=$newpath
   fi
   if [ -n "${DYLD_LIBRARY_PATH}" ]; then
      drop_from_path $DYLD_LIBRARY_PATH ${old_restframessys}/lib
      DYLD_LIBRARY_PATH=$newpath
   fi
   if [ -n "${SHLIB_PATH}" ]; then
      drop_from_path $SHLIB_PATH ${old_restframessys}/lib
      SHLIB_PATH=$newpath
   fi
   if [ -n "${LIBPATH}" ]; then
      drop_from_path $LIBPATH ${old_restframessys}/lib
      LIBPATH=$newpath
   fi
   if [ -n "${PYTHONPATH}" ]; then
      drop_from_path $PYTHONPATH ${old_restframessys}/lib
      PYTHONPATH=$newpath
   fi
   if [ -n "${CPLUS_INCLUDE_PATH}" ]; then
      drop_from_path $CPLUS_INCLUDE_PATH ${old_restframessys}/include
      CPLUS_INCLUDE_PATH=$newpath
   fi
   if [ -n "${CPATH}" ]; then
      drop_from_path $CPATH ${old_restframessys}/include
      CPATH=$newpath
   fi
fi

if [ -z "${PATH}" ]; then
   PATH=$PWD/RestFrames/bin; export PATH
else
   PATH=$PWD/RestFrames/bin:$PATH; export PATH
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
   LD_LIBRARY_PATH=$PWD/RestFrames/lib; export LD_LIBRARY_PATH       # Linux, ELF HP-UX
else
   LD_LIBRARY_PATH=$PWD/RestFrames/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
fi

if [ -z "${DYLD_LIBRARY_PATH}" ]; then
   DYLD_LIBRARY_PATH=$PWD/RestFrames/lib; export DYLD_LIBRARY_PATH   # Mac OS X
else
   DYLD_LIBRARY_PATH=$PWD/RestFrames/lib:$DYLD_LIBRARY_PATH; export DYLD_LIBRARY_PATH
fi

if [ -z "${SHLIB_PATH}" ]; then
   SHLIB_PATH=$PWD/RestFrames/lib; export SHLIB_PATH                 # legacy HP-UX
else
   SHLIB_PATH=$PWD/RestFrames/lib:$SHLIB_PATH; export SHLIB_PATH
fi

if [ -z "${LIBPATH}" ]; then
   LIBPATH=$PWD/RestFrames/lib; export LIBPATH                       # AIX
else
   LIBPATH=$PWD/RestFrames/lib:$LIBPATH; export LIBPATH
fi

if [ -z "${PYTHONPATH}" ]; then
   PYTHONPATH=$PWD/RestFrames/lib; export PYTHONPATH
else
   PYTHONPATH=$PWD/RestFrames/lib:$PYTHONPATH; export PYTHONPATH
fi

if [ -z "${CPLUS_INCLUDE_PATH}" ]; then
   CPLUS_INCLUDE_PATH=$PWD/RestFrames/include; export CPLUS_INCLUDE_PATH
else
   CPLUS_INCLUDE_PATH=$PWD/RestFrames/include:$CPLUS_INCLUDE_PATH; export CPLUS_INCLUDE_PATH
fi

if [ -z "${CPATH}" ]; then
   CPATH=$PWD/RestFrames/include; export CPATH
else
   CPATH=$PWD/RestFrames/include:$CPATH; export CPATH
fi

unset old_restframessys

