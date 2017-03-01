#! /bin/bash
#
# Create the table file in the products directory
#
# Arguments:
#  1 - path to the table file.
#
# This script overwrites previously existing installs of the table file.
#

# Arguments
destination_file=$1

# Discover the versions of the products on which BTrk depends
clhep_ver=`ups active | awk '$1 == "clhep" {print $2}'`
root_ver=`ups active | awk '$1 == "root" {print $2}'`
gcc_ver=`ups active | awk '$1 == "gcc" {print $2}'`

# A table file needs the qualifiers string in three formats:
#   - colon delimited
#   - dot deliminted
#   - colon delimited with plus signs.
# Define them all here and use them in the right spots in the table file.

# qual1:qual2:qual3
# No leading colon on the first qualifer
colon_qualifiers=${COMPILER_CODE}

# qual1.qual2.qual3
# No leading dot on the first qualifer
dot_qualifiers=${COMPILER_CODE}

# +qual1:+qual2:+qual3
# No leading + or colon on the first qualifer
plus_qualifiers=${COMPILER_CODE}

# Write the table file in place
cat > ${destination_file} <<EOF
File    = table
Product = ${PACKAGE_NAME}

#*************************************************
# Starting Group definition
Group:

Flavor     = ANY
Qualifiers = "${colon_qualifiers}:debug"

  Action = GetFQDir
      envSet( \${UPS_PROD_NAME_UC}_FS, "" )
      execute( "get-directory-name subdir", NO_UPS_ENV, \${UPS_PROD_NAME_UC}_FS )
      envSet (BTRK_FQ, \${\${UPS_PROD_NAME_UC}_FS}.${dot_qualifiers}.debug)
      setupRequired( clhep ${clhep_ver} -q +${plus_qualifiers}:+debug )
      setupRequired( root  ${root_ver} -q +${plus_qualifiers}:+debug )

Flavor     = ANY
Qualifiers = "${colon_qualifiers}:prof"

  Action = GetFQDir
      envSet( \${UPS_PROD_NAME_UC}_FS, "" )
      execute( "get-directory-name subdir", NO_UPS_ENV, \${UPS_PROD_NAME_UC}_FS )
      envSet (BTRK_FQ, \${\${UPS_PROD_NAME_UC}_FS}.${dot_qualifiers}.prof)
      setupRequired( clhep ${clhep_ver} -q +${plus_qualifiers}:+prof )
      setupRequired( root  ${root_ver} -q +${plus_qualifiers}:+prof )

Common:
  Action = setup
    prodDir()
    setupEnv()
    envSet (BTRK_INC, \${BTRK_DIR}/include )

    exeActionRequired(GetFQDir)
    envSet (BTRK_LIB, \${BTRK_DIR}/\${BTRK_FQ}/lib )
    envSet (BTRK_VERSION, \${UPS_PROD_VERSION} )

    setupRequired( gcc   ${gcc_ver} )

    if ( test \`uname\` = "Darwin" )
      pathPrepend(DYLD_LIBRARY_PATH, \${\${UPS_PROD_NAME_UC}_LIB})
    else()
      pathPrepend(LD_LIBRARY_PATH, \${\${UPS_PROD_NAME_UC}_LIB})
    endif ( test \`uname\` = "Darwin" )

End:
# End Group definition
#*************************************************
EOF

unset clhep_ver
unset root_ver
unset gcc_ver
