#!/usr/bin/perl
# Runs the Prepartion codes for the grid routing
# 18 May 2012 jroundy@princeton.edu
use strict;
use warnings;
use FileHandle;
use Term::ANSIColor qw(:constants);

my $start = time();
#   Basin related variables
my $basin ="AFR";
my $res = "900s"; #Resolution given in arc seconds

#   Routing related variables
my $vnstr = 2.0; #non stream velocity
my $vstr = 3.0; #stream velocity
my $soid = 1; #cut off for what is a stream (stream order)
my $sid = 1; #1-use slope in calculating velocity for non-streams

#   directory related variables
my $baseDir = "/home/freeze/water_monitor/African_Drought_Monitor/SOURCE/Grid_Routing";
my $GISDir = "$baseDir/GIS/${basin}_$res";

# Model Inputs
my $maskfile = "";
my ($mnx,$mny,$mxmin,$mymin,$mxres,$myres) = (0);
if($basin eq 'AFR')
{
  # Africa
  $maskfile = "$GISDir/AFR_VIC_Mask.bin";
  $mnx = 296;
  $mny = 292;
  $mxmin = -18.875000;
  $mymin = -34.875000;
  $mxres = 0.25000;
  $myres = 0.25000;
}
if($basin eq 'CONUS')
{
  # CONUS Domain
  $maskfile = "/home/wind/jroundy/Workspace/Routing/VIC/mask.bin";
  $mnx = 464;
  $mny = 196;
  $mxmin = -124.9375;
  $mymin = 25.0625;
  $mxres = 0.125000;
  $myres = 0.125000;
}

my $WorkDir = "$baseDir/Workspace";
system("mkdir -p $WorkDir");
chdir($WorkDir);

#Run the Fortran Program
&Prep_Routing($basin,$res);

my $end = time();
print "Time taken was ", ($end - $start), " seconds\n";

###############################################################################
# Runs the fortran program for Preparing the Routing files
sub Prep_Routing
{
    my ($basin,$res) = @_;
    my $prog = "Prep_Routing";
    my $exe = "$prog.exe";
    my $deminfo = `head -6 $GISDir/${basin}_${res}_DIR.ctl | tail -2 | awk '{print \$2,";",\$4,";",\$5}'`;
    my @xy= split("\n",$deminfo);
    my ($nx,$xmin,$xres) = split(";",$xy[0]);
    my ($ny,$ymin,$yres) = split(";",$xy[1]);

    my $fh = new FileHandle (">$prog.nml");
    print $fh <<EOF;
&${prog}_nml
basin='$basin',
outdir='$GISDir/${basin}_${res}',
maskfile='$maskfile',
nx=$nx,
ny=$ny,
xmin=$xmin,
ymin=$ymin,
xres=$xres,
yres=$yres,
mnx=$mnx,
mny=$mny,
vnstr=$vnstr,
vstr=$vstr,
soid=$soid,
sid=$sid,
mxmin=$mxmin,
mymin=$mymin,
mxres=$mxres,
myres=$myres /

EOF
    close $fh;
    my $fname = "$GISDir/${basin}_${res}_Model_Grid.ctl";
    print("$fname\n");
    my $fh2 = new FileHandle (">$fname");
    print $fh2 <<EOF;
dset ${GISDir}/${basin}_${res}_Model_Grid.bin
title Model grid for $basin at $res resolution
options little_endian
undef -999
xdef $nx linear $xmin $xres
ydef $ny linear $ymin $yres
zdef 1 linear 1 1
tdef 1 linear 00Z01Oct1982  1mo
vars 2
ix 0 -1,40,4 X corrdinate of model grid
iy 0 -1,40,4 Y corrdinate of model grid
endvars
EOF
    close $fh2;

    $fname = "$GISDir/${basin}_${res}_ADIR.ctl";

    $fh2 = new FileHandle (">$fname");
    print $fh2 <<EOF;
dset ${GISDir}/${basin}_${res}_ADIR.bin
title x,y direction for $basin at $res resolution adjusted for routing
options little_endian
undef -9999
xdef $nx linear $xmin $xres
ydef $ny linear $ymin $yres
zdef 1 linear 1 1
tdef 1 linear 00Z01Oct1982  1mo
vars 2
xdir 0 -1,40,4 X corrdinate of flow for gridcell
ydir 0 -1,40,4 Y corrdinate of flow for gridcell
endvars
EOF
    close $fh2;

    $fname = "$GISDir/${basin}_${res}_VEL.ctl";

    $fh2 = new FileHandle (">$fname");
    print $fh2 <<EOF;
dset ${GISDir}/${basin}_${res}_VEL.bin
title grid velocity $basin at $res resolution
options little_endian
undef -999.9
xdef $nx linear $xmin $xres
ydef $ny linear $ymin $yres
zdef 1 linear 1 1
tdef 1 linear 00Z01Oct1982  1mo
vars 1
vel 0 99 velocity parameter [m/s]
endvars
EOF
    close $fh2;

    $fname = "$GISDir/${basin}_${res}_ARE.ctl";

    $fh2 = new FileHandle (">$fname");
    print $fh2 <<EOF;
dset ${GISDir}/${basin}_${res}_ARE.bin
title Area of grids for $basin at $res resolution
options little_endian
undef -999.9
xdef $nx linear $xmin $xres
ydef $ny linear $ymin $yres
zdef 1 linear 1 1
tdef 1 linear 00Z01Oct1982  1mo
vars 1
area 0 99 grid area [m2]
endvars
EOF
    close $fh2;

    $fname = "$GISDir/${basin}_${res}_LEN.ctl";

    $fh2 = new FileHandle (">$fname");
    print $fh2 <<EOF;
dset ${GISDir}/${basin}_${res}_LEN.bin
title Distance to next grid for $basin at $res resolution
options little_endian
undef -999.9
xdef $nx linear $xmin $xres
ydef $ny linear $ymin $yres
zdef 1 linear 1 1
tdef 1 linear 00Z01Oct1982  1mo
vars 1
len 0 99 length [m]
endvars
EOF
    close $fh2;

    # Run Code
    my $logfile = "log_GIS_Process.$$";
    $logfile = "$WorkDir/$logfile";
    system("$baseDir/src/Prep_Routing/$exe");
#     system("$baseDir/src/Prep_Routing/$exe >> $logfile 2>&1  ");
}
