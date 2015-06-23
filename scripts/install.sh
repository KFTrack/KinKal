#! /bin/bash
#
# Copy the files from the source and build areas to create a UPS product.
#
# This script deletes previously existing installs of the same product+version.
# Use with care.
#

# Check that the installation directoy has been defined.
if [ "${PRODUCTS_INSTALL}" = '' ];then
    echo "The environment variable PRODUCTS_INSTALL is not set."
    echo "You must define where to install the products before sourcing this script."
    return 1
fi

# Build the names of the directories into which we will write things
flavour=`ups flavor`
fq=${flavour}-${COMPILER_CODE}-${DEBUG_LEVEL}
topdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}
proddir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}
verdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}.version
libdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/${fq}
incdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/include
upsdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/ups

# I am not sure what this file is properly called.
# I am calling it the fqfile, which is short for flavor qualifier file.
fqfile=${verdir}/${flavour}_${COMPILER_CODE}_${DEBUG_LEVEL}

# Make directories, if needed.
if ! [ -e ${topdir} ];then
  mkdir ${topdir}
fi

if ! [ -e ${proddir} ];then
  mkdir ${proddir}
fi

if ! [ -e ${verdir} ];then
  mkdir ${verdir}
fi

# Remove all content that we will recreate.
# In all cases create required directories.
if ! [ -e ${libdir} ];then
  mkdir ${libdir}
else
  /bin/rm -rf ${libdir}
  mkdir ${libdir}
fi

if ! [ -e ${incdir} ];then
  mkdir ${incdir}
else
  /bin/rm -rf ${incdir}
  mkdir ${incdir}
fi

if ! [ -e ${upsdir} ];then
  mkdir ${upsdir}
else
  /bin/rm -rf ${upsdir}
  mkdir ${upsdir}
fi

if [ -e ${fqfile} ];then
  /bin/rm -rf ${fqfile}
fi

# FIXME: change tar to rsync to preserve times?
# Copy the required parts of the source directory to the installation area.
tar cf - --exclude-from  ${PACKAGE_SOURCE}/etc/tarExcludePatterns.txt \
    -C ${PACKAGE_SOURCE} ${PACKAGE_NAME} | tar xf - -C ${proddir}/include

tar cf - --exclude-from  ${PACKAGE_SOURCE}/etc/tarExcludePatterns.txt \
    -C ${PACKAGE_SOURCE} ups | tar xf - -C ${proddir}

# Copy the lib directory to the installation area
tar cf - lib  | tar xf - -C ${libdir}

# Fixme: temporary hack until #include directives are cleaned up.
ln -s ${incdir}/BTrk ${incdir}/BaBar

# Create the ups fq file.
cat > ${fqfile} <<EOF
FILE = version
PRODUCT = ${PACKAGE_NAME}
VERSION = ${PACKAGE_VERSION}

#*************************************************
#
FLAVOR = ${flavour}
QUALIFIERS = "${COMPILER_CODE}:${DEBUG_LEVEL}"
  DECLARER = `whoami`
  DECLARED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  MODIFIER = `whoami`
  MODIFIED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  PROD_DIR = ${PACKAGE_NAME}/${PACKAGE_VERSION}
  UPS_DIR = ups
  TABLE_FILE = ${PACKAGE_NAME}.table

EOF
