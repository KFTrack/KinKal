# Refactoring BField Interface

Create your own ups products directory:

1. Create an empty directory. No part of the name is siginficant but it seems sensible
   to adopt the convention that the last part of the name is artexternals.

2. cd to the directory you just made

3. cp -r  /cvmfs/mu2e.opensciencegrid.org/artexternals/.upsfiles .
   You can copy this from any existing ups repository

4. On every login, after you do "setup mu2e"  you need to define the following environment variables
   export PRODUCTS_INSTALL=/absolute/path/to/your/products/directory
   export PRODUCTS=${PRODUCTS}:${PRODUCTS_INSTALL}

We can work on structuring step 4 differently - maybe if PRODUCTS_INSTALL is defined then have "setup mu2e" do
the addition to PRODUCTS?

-----------------------------------------------------------

Setting up the source directory:

1. setup mu2e
2. Step 4 from creating your own products directory
3. cd to a clean directory
4. git clone https://github.com/KFTrack/BTrk.git
    - This will make a subdirectory BTrk
5. If you wish to build the head of BTrk, proceed to the
   next section. If you wish to build a tag of BTrk, then
      git checkout -b work tag_name
   Here work is an arbitrary name for your working branch
   and tag_name is the name of the tag you wish to build.

Making a prof build
1. cd to a clean directory that must be outside of the source directory tree.
   ( This examples assumes that you use same directory that holds BTrk ).
2. mkdir build_prof
3. cd build_prof
4. source ../BTrk/scripts/newBuild.sh prof
   - this will but two files in the directory: SConstruct, setup.sh
   - the prof option is now hardcoded for all future builds.
5. source setup.sh
6. scons -j 8
    - A recent build took 1 min 30 sec of real time and 5 minutes of CPU time.
    - this no longer writes anything in the source tree
    - .os files are written to a new subdirectory, BTrk/
    - libraries are written to a new subdirectory, lib/
7. source ../BTrk/scripts/install.sh
    - this will write the headers and libraries to the products directory made above

As you develop code, edit in the source directory and build in the build directory.
the build scripts will never write into the source directory.

On subsequent logins if you want to rebulid:

1. setup mu2e
2. Step 4 from creating your own products directory
3. cd build_prof
4. source setup.sh

Making a debug build

1. If you have just done the prof build it's easiest to log out and log in again.
2. setup mu2e
3. Step 4 from creating your own products directory
4. cd to the a clean directory - same guidance as above
5. mkdir build_debug
6. cd build_debug
7. source ../BTrk/scripts/newBuild.sh debug
8. Repeat steps 5, 6, 7 from building the prof build.

Some ideas:
1. You may wish to put your source directory on a backed up disk, such as /grid/fermiapp/mu2e or your home disk,
   put your build directories on a scratch disk and put the products are on /mu2e/app, which is not backed up but which
   is mounted on the GPGrid worker nodes.
2. To run offsite We will need user products areas on cvmfs.


Building Mu2e Offline:

1. setup mu2e
2. Step 4 from creating your own products directory
    - You only need to update the definition of PRODUCTS.
    - The definition of PRODUCT_INSTALL is not needed
3. git clone ssh://p-mu2eofflinesoftwaremu2eoffline@cdcvs.fnal.gov/cvs/projects/mu2eofflinesoftwaremu2eoffline/Offline.git
4. scons -j 8

Running art:

1. Add the following line to the services block of any fcl file that you use:
     BTrkHelper : @local::BTrkHelperDefault

What did I change:

In BTrk
1. created subdirectory BTrk and moved all source directories down one level
2. Added SConstruct, setup.sh, and directories: etc/, python/, scripts/, ups/
3. The big change is that I refactored SConstruct. It now is configured to do
   an out-of-source build. THe helper functions are in python//helpers.py
4. Edited all SConcript files
    - the helper classname is changed to build_helper from mu2e_helper
    - I removed the module and data product builds since they are not used.

In Offline:
1. setup.sh
    - Use the BTrk ups product
    - remove section about BaBar svn code
2. SConstruct
    - Use the BTrk ups product
3. Remove the BaBar directory
4. BTrkHelper/src/BTrkHelper_service.cc
     - fixed typo
5. Put copies of material description files into BTrkHelper/data
     ElementsList.data
     IsotopesList.data
     MaterialsList.data
   Edit BTrkHelper/fcl/prolog.fcl to use these files
6. Analyses/test/genReco.fcl
     - add configuration for the BTrkHelper service

When you wish to make a new version of BTrk, without changing any of
the products on which it depends:

 1. In setup.sh, edit PACKAGE_VERSION to give it the new value
 2. commmit changes, tag, push changes and tag
 3. build and install both debug and prof versions
 4. If necessary, edit Offline/setup.sh to use the new version


If you update to a new compiler, or a new version/qualifier of any of the
products that BTrk sets up, you need to:
1. Edit setup.sh
    a) Change the setup lines as needed; see below for details on
       upgrading to match a new version of art.
    b) If you change the compiler, edit the COMPILER_CODE.
    c) You should also bump PACKAGE_VERSION.  For now, if you want to
       change any UPS dependencies, you need a new version of BTrk. We
       have not yet establish a policy to maintain this bookkeeping
       using UPS qualifiers.
2. Edit installTableFile.sh
    a) The following are dealt with automatically by installTableFile.sh
          - the version of each ups product
          - the compiler code qualfiier, which is common to all products
    b) If you add other qualifers in step 1., update the variables
       colon_qualifiers, dot_qualifiers, plus_qualifiers

---------------------------------------------------------------------------
Upgrading for compatibility with a new version of art:

0.
These instructions guide you through step 1. in the previous section
for the example of upgrading to a new version of art.

1.
Setup the new version of art. For example:
 setup art v1_17_02 -qe9:prof

2.
List all of the products that have been setup by art:
ups active

3.
Edit setup.sh and find all ups setup commands; change the versions
and qualifiers to match those shown by ups active.

4.
In step 1., in the qualifiers to setup art, there is a "e" qualifier,
"e9" in the above example. In setup.sh, edit the definition of
COMPILER_CODE to match.

5.
Discover which versions of scons are available:
ups list -aK+ scons
This is not set up by setting up art so you will not see it in the
listing generated by ups active. In general choose the most recent version.
Edit setup.sh to choose this version.

6.
The policy set by the art team is that the ups product cetpkgsupport has
a version declared "current".  Our policy is to use this one rather than
to choose specific versions - there are no binary compatibility issues with
this product.  So it is not necessary to edit setup.sh.

7.
In setup.sh, edit the definition of PACKAGE_VERSION to match the tag
you have just made (or will make as soon as tests pass).

8.
If there were any changes in qualifiers, follow the instructions
in the previous section.

### Quick Rundown of commands for building local BTrk and integrating with Offline:
1. build_prof (compiling BTrk changes)
  * `source setup.sh`
  * `scons -j 8`
  * `source ../BTrk/scripts/install.sh`
2. Offline (for compiling offline changes)
  * `cd Offline`
  * `source setup.sh`
  * `scons -j 8`
3. test (for running validation script via TrkAna)
  * `cd test`
  * `source ../Offline/setup.sh`
  * `mu2e --config TrkDiag/test/TrkAna.fcl --source sim.owner.cd3-beam-g4s4-detconversion.version.sequencer.art --TFile MYOUTPUT.root --nevts=100 >| MYLOG.log`
