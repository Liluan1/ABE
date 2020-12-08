Dmitri here are the main things left to do:

1) remove third-party (I don't think we need it?)

2) rework CMake to that we use the  installed palsade library. 
use the CMakeLists.User.txt approach. I had to hardwire the library in the cmake files because I did not understand the User.txt file was there. 


3) this is called abe but it has abe signatures and trapdoor. We want to split this out into three separate repos called abe signatures and trapdoor, each with their own src/abe. src/signatures and src/trapdoor

4) update the Readme.md file with new build instructions. palisade 1.10.6 should be installed on the system first, then this code should be made.

Dave
