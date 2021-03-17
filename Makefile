LFLAGS = -L$(ROOTSYS)/lib -lCore -lRIO -lEG -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui -lMinuit -lm -ldl

TARGET = corsikaConverter

default: $(TARGET)

$(TARGET): $(TARGET).cc
	g++ -O2 -Wno-deprecated -std=c++17 $(TARGET).cc -o corsikaConverter -pthread -I$(ROOTSYS)/include $(LFLAGS) -pthread -rdynamic -Wall


clean:
	$(RM) corsikaConverter
