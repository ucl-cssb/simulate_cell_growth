SRC_DIR = src
BIN_DIR = bin


sim:
	@mkdir -p ${BIN_DIR}
	g++ -std=gnu++11 $(SRC_DIR)/cell.cpp $(SRC_DIR)/sample.cpp $(SRC_DIR)/main.cpp -o $(BIN_DIR)/sim

clean:
	rm $(BIN_DIR)/sim
