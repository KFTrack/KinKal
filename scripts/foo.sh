#
# Original author Rob Kutschke
#
# Setup the environment to build standalone BaBar code.
#
# The one required argument is the debug level, either prof or debug.
# It is sticky.  To build both debug and prof you need to create
# two build areas and run this script once in each area.
#

if [ "`basename $0 2>/dev/null`" = "setup.sh" ];then
    echo "You should be sourcing this file, not executing it."
    exit 1
fi

# The root of the source to be built is the directory in which this script is found.
package_source=`cd "$(dirname ${BASH_SOURCE})" >/dev/null 2>&1 && /bin/pwd | sed 's/\/scripts//' `
echo $package_source
basename $package_source

