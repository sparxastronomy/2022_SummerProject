# Script to obtain catalog information from the EagleDB at z=0
echo "Getting StrongFBL0025N0376_FoF"
wget --http-user=pnm600 --http-passwd=CK205HQn --cookies=on --keep-session-cookies --save-cookies=cookie.txt --load-cookies=cookie.txt -O StrongFBL0025N0376_FoF.csv "http://galaxy-catalogue.dur.ac.uk:8080/Eagle?action=doQuery&SQL=select GroupCentreOfPotential_x as X, GroupCentreOfPotential_y as Y, GroupCentreOfPotential_z as Z, GroupMass as Mass from Physics_vars..StrongFBL0025N0376_FoF where SnapNum=28"
echo "Getting WeakFBL0025N0376_FoF"
wget --http-user=pnm600 --http-passwd=CK205HQn --cookies=on --keep-session-cookies --save-cookies=cookie.txt --load-cookies=cookie.txt -O WeakFBL0025N0376_FoF.csv "http://galaxy-catalogue.dur.ac.uk:8080/Eagle?action=doQuery&SQL=select GroupCentreOfPotential_x as X, GroupCentreOfPotential_y as Y, GroupCentreOfPotential_z as Z, GroupMass as Mass from Physics_vars..WeakFBL0025N0376_FoF where SnapNum=28"
echo "Getting ViscHiL0050N0752_FoF"
wget --http-user=pnm600 --http-passwd=CK205HQn --cookies=on --keep-session-cookies --save-cookies=cookie.txt --load-cookies=cookie.txt -O ViscHiL0050N0752_FoF.csv "http://galaxy-catalogue.dur.ac.uk:8080/Eagle?action=doQuery&SQL=select GroupCentreOfPotential_x as X, GroupCentreOfPotential_y as Y, GroupCentreOfPotential_z as Z, GroupMass as Mass from Physics_vars..ViscHiL0050N0752_FoF where SnapNum=28"
echo "Getting ViscLoL0050N0752_FoF"
wget --http-user=pnm600 --http-passwd=CK205HQn --cookies=on --keep-session-cookies --save-cookies=cookie.txt --load-cookies=cookie.txt -O ViscLoL0050N0752_FoF.csv "http://galaxy-catalogue.dur.ac.uk:8080/Eagle?action=doQuery&SQL=select GroupCentreOfPotential_x as X, GroupCentreOfPotential_y as Y, GroupCentreOfPotential_z as Z, GroupMass as Mass from Physics_vars..ViscLoL0050N0752_FoF where SnapNum=28"
echo "Getting RefL0100N1504_FOF"
wget --http-user=pnm600 --http-passwd=CK205HQn --cookies=on --keep-session-cookies --save-cookies=cookie.txt --load-cookies=cookie.txt -O RefL0100N1504_FOF.csv "http://galaxy-catalogue.dur.ac.uk:8080/Eagle?action=doQuery&SQL=select GroupCentreOfPotential_x as X, GroupCentreOfPotential_y as Y, GroupCentreOfPotential_z as Z, GroupMass as Mass from Eagle..RefL0100N1504_FOF where SnapNum=28"
