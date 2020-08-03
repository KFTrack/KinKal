#
# Original author Rob Kutschke
#
# Setup the environment to build standalone BaBar code.
#
# There is one required argument, the debug level, which must be
# one of prof or debug.
#
# Prior to calling this the user must set the environment variable
# PACKAGE_SOURCE with the path to the checked-out source of the
# package to be built.
#

#if [ "`basename $0 2>/dev/null`" = "setup.sh" ];then
#    echo "You should be sourcing this file, not executing it."
#    return 1
#fi

if [ "${PACKAGE_SOURCE}" = '' ];then
    echo "The environment variable PACKAGE_SOURCE is not set or is an empty string."
    echo "You must set this environment variable to use this script."
    return 1
fi

if ! [ -d ${PACKAGE_SOURCE} ];then
   echo "The directory named by the environment variable PACKAGE_SOURCE does not exist."
   echo "PACKAGE_SOURCE: " ${PACKAGE_SOURCE}
   return 1
fi

# The target directory for the build is the directory from which the
# source command is issued.
export BUILD_BASE=`pwd -P`
echo "Source to be built: " ${PACKAGE_SOURCE}
echo "Root of build:      " ${BUILD_BASE}
echo "Debug level:        " ${DEBUG_LEVEL}

# Tell SConstruct where to find helpers.py
if [ "${PYTHONPATH}" = '' ];then
 export PYTHONPATH=${PACKAGE_SOURCE}/python
else
 export PYTHONPATH=${PYTHONPATH}:${PACKAGE_SOURCE}/python
fi

# Tell python not to write out byte code files.
# Fixme: Otherwise they write byte code files into the source directory.
# I would prefer to specify a cache inside the build directory but
# I don't know how to do that.
export PYTHONDONTWRITEBYTECODE=1

export ROOT_INC=$(root-config --incdir) #${ROOTSYS}/include
export ROOT_LIB=$(root-config --libdir) #${ROOTSYS}/lib

export GCC_VERSION=`gcc -dumpversion | sed 's/\./_/g'`
