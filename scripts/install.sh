#! /bin/bash
#
# Copy the files from the source and build areas to create a UPS product.
#
# This script unconditionally deletes previously existing installs of the
# same product+version+qualifiers: use with care.
#

# Check that the installation directoy has been defined.
if [ "${PRODUCTS_INSTALL}" = '' ];then
    echo "The environment variable PRODUCTS_INSTALL is not set."
    echo "You must define where to install the products before sourcing this script."
    return 1
fi

# There are two representations of operating system UPS flavor:
# old style, for example: Linux64bit+2.6-2.12_e7
# new style, for example: slf6.x86_64
# We need them both.
old_flavour=`ups flavor`
new_flavour=`get-directory-name subdir`

# Build the names of the directories into which we will write things
fq=${new_flavour}.${COMPILER_CODE}.${DEBUG_LEVEL}
topdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}
proddir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}
verdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}.version
libdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/${fq}
incdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/include
upsdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/ups

# I am not sure what this file is properly called.
# I am calling it the fqfile, which is short for flavor qualifier file.
fqfile=${verdir}/${old_flavour}_${COMPILER_CODE}_${DEBUG_LEVEL}

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
# (Re)create required directories.
if [ -e ${libdir} ];then
  /bin/rm -rf ${libdir}
fi
mkdir ${libdir}

if [ -e ${incdir} ];then
  /bin/rm -rf ${incdir}
fi
mkdir ${incdir}

if [ -e ${upsdir} ];then
  /bin/rm -rf ${upsdir}
fi
mkdir ${upsdir}

if [ -e ${fqfile} ];then
  /bin/rm -rf ${fqfile}
fi

# Copy the required parts of the source directory to the installation area:

# Header files:
rsync -ar --exclude-from  ${PACKAGE_SOURCE}/etc/tar_exclude_for_include.txt \
    ${PACKAGE_SOURCE}/${PACKAGE_NAME} ${proddir}/include

# UPS table file
rsync -ar --exclude-from  ${PACKAGE_SOURCE}/etc/tar_exclude_for_include.txt \
    ${PACKAGE_SOURCE}/ups ${proddir}

# Libaries
rsync -ar lib ${libdir}

# A copy of the full source
rsync -ar --exclude-from  ${PACKAGE_SOURCE}/etc/tar_exclude_for_source.txt \
  ${PACKAGE_SOURCE}/${PACKAGE_NAME} ${proddir}/source

# Create the ups fq file.
cat > ${fqfile} <<EOF
FILE = version
PRODUCT = ${PACKAGE_NAME}
VERSION = ${PACKAGE_VERSION}

#*************************************************
#
FLAVOR = ${old_flavour}
QUALIFIERS = "${COMPILER_CODE}:${DEBUG_LEVEL}"
  DECLARER = `whoami`
  DECLARED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  MODIFIER = `whoami`
  MODIFIED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  PROD_DIR = ${PACKAGE_NAME}/${PACKAGE_VERSION}
  UPS_DIR = ups
  TABLE_FILE = ${PACKAGE_NAME}.table

EOF

unset old_flavour
unset new_flavour
unset fq
unset topdir
unset proddir
unset verdir
unset libdir
unset incdir
unset upsdir
unset fqfile
