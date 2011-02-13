Welcome to the second release of the Data Collaboration toolbox. 

Throughout this readme, ./DClabV2/ indicates the base directory of
your installation (e.g., C:/program files/DClabV2/).

Installation instructions are available in html format by changing
into the ./DClabV2/docs/html/ directory and launching manual.html. The
impatient user may follow the instructions given below.


                     ==To Install==

Open MATLAB, and make ./DClabV2/ your current directory. At the
command prompt type 'DCsetup'. This will add a handful of directories
to your path. To make the installation permanent, use file | set path
in the MATLAB file menu, or manually add these directories (viewable
by opening DCsetup.m) to your pathdef.m or startup.m file. You may
then test the installation by typing 'dctest'. If you are working in
unix, you may be prompted to install the SeDuMi binaries. Hopefully
this will consist of simply executing
'./DClabV2/SeDuMi_1_1/install_sedumi'. This script mexes
several .c files and consequently its sucessful execution will depend
on your system's configuration. If it fails to execute properly, you
use unix, you know the drill, fix it.

                  ===What to do next===

The best place to start is the html documentation, available by
launching ./DClabV1/docs/html/manual.html in a browser window.

After that, you're pretty much on your own. The directory
./DClabV2/fileTemplates/ contains a few templates to get you
started.

For questions and especially suggestions, email Trent Russi at
trussi@berkeley.edu

Good Luck!
