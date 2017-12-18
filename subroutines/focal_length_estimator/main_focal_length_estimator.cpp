//
// Created by danielbord on 12/18/17.
//

int main(int argc, char *argv[]) {

    namespace po = boost::program_options;
    try {

        po::options_description desc("Optimization input options");
        desc.add_options()
                ("help", "Print help message")
                ("fund_f", po::value<std::string>(&f_f)->required(), "File with %n_pic estimated fundamental matrices");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);


        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return -1;
        }
        boost::program_options::notify(vm);

    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -2;
    }


    return 0;
}