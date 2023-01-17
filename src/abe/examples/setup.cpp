#include "abecontext.h"
#include "utils/serial.h"

using namespace lbcrypto;

// const std::string DATAFOLDER = "serial-data";

void parse_args(int argc, char *argv[], usint &numAttributes,
                std::string &mpkFile, std::string &mskFile) {
  //   if (argc < 4) {
  //     std::cout << "Usage: " << argv[0] << " <numAttributes>"
  //               << " <mpkFile>"
  //               << " <mskFile>" << std::endl;
  //     exit(1);
  //   }
  //   numAttributes = atoi(argv[1]);
  //   mpkFile = argv[2];
  //   mskFile = argv[3];
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {
      numAttributes = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-mpk") == 0) {
      mpkFile = argv[++i];
    } else if (strcmp(argv[i], "-msk") == 0) {
      mskFile = argv[++i];
    }
  }
}

int main(int argc, char *argv[]) {
  usint ringsize = 1024;
  usint numAttributes;
  usint base = 64;
  std::string mpkFile;
  std::string mskFile;
  parse_args(argc, argv, numAttributes, mpkFile, mskFile);

  std::cout << "Used parameters:" << std::endl;
  std::cout << "Ring size: " << ringsize << std::endl;
  std::cout << "Number of attributes: " << numAttributes << std::endl;
  std::cout << "Base: " << base << std::endl << std::endl;

  ABEContext<NativePoly> context;
  std::cout << "Generating a context under these parameters" << std::endl
            << std::endl;
  context.GenerateCPABEContext(numAttributes, ringsize, base);

  std::cout << "Generating master secret key and master public key"
            << std::endl;

  // Generate master keys

  CPABEMasterPublicKey<NativePoly> mpk;
  CPABEMasterSecretKey<NativePoly> msk;
  context.Setup(&mpk, &msk);

  //   std::cout << "Setup: " << duration << " ms" << std::endl << std::endl;

  // Serialize master public keys
  if (!Serial::SerializeToFile(mpkFile, mpk, SerType::BINARY)) {
    std::cerr << "Error writing serialization of the master public key to "
              << mpkFile << std::endl;
    return 1;
  }

  // Serialize master secret keys
  if (!Serial::SerializeToFile(mskFile, msk, SerType::BINARY)) {
    std::cerr << "Error writing serialization of the master secret key to "
              << mskFile << std::endl;
    return 1;
  }
  return 0;
}