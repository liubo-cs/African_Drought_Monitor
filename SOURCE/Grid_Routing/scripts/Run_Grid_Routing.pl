#!/usr/bin/perl
# Runs the Grid Routing 
# Original - 18 May 2012 jroundy@princeton.edu
# Updated - 21 Aug 2013 jroundy@princeton.edu
use strict;
use warnings;
use FileHandle;

my $start = time();

# Domain Inputs
my $basin ="AFR";
my $res = "900s"; #Resolution given in arc seconds

# Directory Inputs
#my $baseDir = "/home/raid20/jroundy/Grid_Routing";
my $baseDir = "/home/freeze/water_monitor/African_Drought_Monitor/SOURCE/Grid_Routing";
my $GISDir = "$baseDir/GIS/${basin}_$res";

# Time Inputs
my $syr = $ARGV[0];#1948; #Start year
my $smn = $ARGV[1];#1; #Start month of simulation
my $sdy = $ARGV[2];#1; #Start Day of simulation
my $nt = $ARGV[3];#24; #Number of time steps to run (see tid)
my $tid = 0;#1; # units of nt, 0-days, 1-months
my $ist = $ARGV[4];#-1; #1-read in previous state, -1 start with zero
my $ost = 1; #output state file 0-none,1-everyday,2-last day of month,3-end of run

# Algorithm Inputs
my $nth = 8; #number of threads to run
my $mxg = 20; #maxium number of cells that the flow can be routed through, depends on resolution of topography and velocity  (300-30s),(20-900s)

# Model Inputs
my $modDir = $ARGV[5];
my ($mtd,$nvr,$rid,$bid,$mnx,$mny,$mxmin,$mymin,$mxres,$myres) = (0);
if($basin eq 'AFR')
{
  # Africa
  #$modDir = "/home/raid20/nchaney/ADM_FOLLOW/Africa_Drought_Monitor/Data/Historical/VIC_output/output_grid_";
  #$modDir = "/home/freeze/water_monitor/African_Drought_Monitor/DATA/VIC/OUTPUT/DAILY/output_grid_";
  $mtd = 2; #Model output time step, 1-monthly, 2-daily
  $nvr = 21; # number of variables in the model output file
  $rid = 3; # variable number for runoff
  $bid = 4; # variable number for baseflow
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
  $modDir = "/home/wind/jroundy/Workspace/Routing/VIC/output_grid_";
  $mtd = 1; #Model output time step, 1-monthly, 2-daily
  $nvr = 21; # number of variables in the model output file
  $rid = 3; # variable number for runoff
  $bid = 4; # variable number for baseflow
  $mnx = 464;
  $mny = 196;
  $mxmin = -124.9375;
  $mymin = 25.0625;
  $mxres = 0.125000;
  $myres = 0.125000;
}

my $WorkDir = "$baseDir/Workspace";
my $streamDir = $ARGV[6];#"$WorkDir/${basin}_$res/Flow";
my $stateDir = $ARGV[7];#"$WorkDir/${basin}_$res/State";

system("mkdir -p $streamDir");
system("mkdir -p $stateDir");
chdir($WorkDir);

#Run the Fortran Program
&Grid_Routing($basin,$res);

my $end = time();
print "Time taken was ", ($end - $start), " seconds\n";

###############################################################################
#
#   Subroutine to Run the Grid_Routing Program
sub Grid_Routing
{
    my ($basin,$res) = @_;
    my $prog = "Grid_Routing";
    my $exe = "$prog.exe";
    my $deminfo = `head -6 $GISDir/${basin}_${res}_DIR.ctl | tail -2 | awk '{print \$2,";",\$4,";",\$5}'`;
    my @xy= split("\n",$deminfo);
    my ($nx,$xmin,$xres) = split(";",$xy[0]);
    my ($ny,$ymin,$yres) = split(";",$xy[1]);

    my $fh = new FileHandle (">$prog.nml");
    print $fh <<EOF;
&${prog}_nml
basin='${basin}_$res',
streamdir='$streamDir',
statedir='$stateDir',
indir='$GISDir/${basin}_${res}',
modDir='$modDir',
nx=$nx,
ny=$ny,
xmin=$xmin,
ymin=$ymin,
xres=$xres,
yres=$yres,
mnx=$mnx,
mny=$mny,
mxmin=$mxmin,
mymin=$mymin,
mxres=$mxres,
myres=$myres,
syr=$syr,
smn=$smn,
sdy=$sdy,
nt=$nt,
tid=$tid,
ist=$ist,
ost=$ost,
mtd=$mtd,
nvr=$nvr,
rid=$rid,
bid=$bid,
nth=$nth,
mxg=$mxg /

EOF
    close $fh;

#     my $fname = "$OutDir/${basin}_${res}_Streamflow.ctl";
#     print("$fname\n");
#     my $fh2 = new FileHandle (">$fname");
#     print $fh2 <<EOF;
# dset ${OutDir}/${basin}_${res}_Streamflow_%y4%m2%d2.bin
# title x,y direction for $basin at $res resolution
# options template little_endian
# undef -999.9
# xdef $nx linear $xmin $xres
# ydef $ny linear $ymin $yres
# zdef 1 linear 1 1
# tdef 31 linear 00Z01JAN1948  1dy
# vars 1
# str 0 99 Routed streamflow
# endvars
# EOF
#     close $fh2;

    # Run Code
    my $logfile = "log_GIS_Process.$$";
    $logfile = "$WorkDir/$logfile";
    system("$baseDir/src/Grid_Routing/$exe ");
#     system("$baseDir/src/Grid_Routing/$exe >> $logfile 2>&1  ");
}
