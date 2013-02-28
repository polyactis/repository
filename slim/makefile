# 2013.2.28 makefile added by Yu Huang (postdoc at ucla)

c++tg=slim

CFLAGS=-O3 -Wl,--no-as-needed 
# -Wl, --no-as-needed  is suggested as a workaround in GSL linking error
# http://ubuntuforums.org/archive/index.php/t-1858271.html
# http://gnu-autoconf.7623.n7.nabble.com/Autotools-g-linking-errors-and-the-no-as-needed-option-td10428.html
LibFlags=-lgsl -lgslcblas -lm

all:	$(c++tg)


.c.o:
	-gcc $(CFLAGS) -c $< -o $@

.cpp.o:
	-g++ $(CFLAGS) -c $< -o $@

$(c++tg): 	%:	%.o
	-g++ $(CFLAGS) $(LibFlags) $(IDIR) $< -o $@ 

$(ctg):	%:	%.o
	-gcc $(CFLAGS) $(LibFlags)  $< -o $@




clean:
	-rm -f *~  *.bb *.bbg *.da *.gcov *.o $(ctg) $(c++tg)
