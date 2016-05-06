/**
 * @file slom/tools/MakeOptions.hpp
 * @brief Brief description
 * 
 */


#ifndef MAKEOPTIONS_HPP_
#define MAKEOPTIONS_HPP_


#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>



#include "../../mtk/build_manifold.hpp"



namespace SLOM {


template<class X>
inline X convert(const std::string& str) {
	X ret;
	std::istringstream oss(str);
	oss >> ret;
	return ret;
}
template<class X>
inline const X& convert(const X&x) {
	return x;
}
inline const std::string& convert_str(const std::string& str) {
	return str;
}

template<class X>
inline std::string convert_str(const X&x) {
	std::ostringstream oss;
	oss << x;
	return oss.str();
}


namespace po=boost::program_options;


struct OptionContainer {
	
	boost::shared_ptr<po::options_description> boost_options;
	
	virtual ~OptionContainer() {}
	
	OptionContainer(const std::string &caption) : boost_options(new po::options_description(caption)) {}
	
	/**
	 * Method which registers members to the options_description.
	 * Must be implemented by the user.
	 */
//	virtual void init() = 0; 
	
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
};


struct MainOptions {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	po::options_description cmd_line, cfg_file;
	
	MainOptions()  {
		cmd_line.add_options()
      ("version,v", "print version string")
      ("help", "produce help message");
	}
	po::variables_map vm;
	
	int parseCommandLine(int argc, char** argv) {
		try{
			
			po::store(po::parse_command_line(argc, argv, cmd_line), vm);
			po::notify(vm);
			
		} catch (const std::exception &e) {
			std::cerr << e.what() << std::endl;
			return -1;
		}
		if (vm.count("help")) {
			std::cout << cmd_line << "\n";
			return 1;
		}
		if (vm.count("version")) {
			std::cout << "Version 0.0.0" << "\n"; // FIXME version string
			return 1;
		}
		return 0;
	}
	
	int parseConfigStream(std::istream& stream) {
		if(!stream.good()) {
			return -2;
		}
		try {
			po::store(po::parse_config_file(stream, cfg_file, true), vm);
			po::notify(vm);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			return -1;
		}
		return 0;
	}
	
	int parseConfigFile(const std::string& filename) {
		std::ifstream file(filename.c_str());
		return parseConfigStream(file);
	}
	
	void storeConfigStream(std::ostream& /*stream*/) {
		// TODO check if that works:
		// http://stackoverflow.com/questions/4703134/is-there-a-way-to-print-config-file-for-boost-program-options
		assert(false && "Not implemented!");
//		boost::property_tree::write_ini(stream, vm);
	}
	
	void storeConfigFile(const std::string& filename) {
		std::ofstream file(filename.c_str());
		storeConfigStream(file);
	}
	
	MainOptions& add(const po::options_description& desc) {
		cmd_line.add(desc);
		cfg_file.add(desc);
		return *this;
	}
	
	MainOptions& add(const OptionContainer& oc) {
		return add(*oc.boost_options);
	}
};

}  // namespace SLOM


#define SLOM_MAKE_OPTIONS(name, members) \
struct name : public SLOM::OptionContainer { \
	typedef name self_; \
	MTK_TRANSFORM(DECLARE_OPTION, members) \
	name() : OptionContainer(#name) { \
		init(); \
	} \
private: \
	void init() { \
		boost_options->add_options() \
		MTK_TRANSFORM(REGISTER_OPTION, members) \
		; \
		std::istringstream dummy; SLOM::po::variables_map vm;\
		SLOM::po::store(SLOM::po::parse_config_file(dummy, *boost_options), vm); SLOM::po::notify(vm);\
	} \
};


#define SLOM_COLLECT_OPTIONS(name, members) SLOM_EXTEND_OPTIONS(SLOM::MainOptions, name, members)

#define SLOM_EXTEND_OPTIONS(base, name, members) \
struct name : public base { \
	MTK_TRANSFORM(ADD_OPTION, members) \
	name() { \
		MTK_TRANSFORM(INSERT_OPTION, members) \
	} \
};


#define ADD_OPTION(type_, name) \
		type_ name;
#define INSERT_OPTION(type, name) add(name);

/**
 * @internal
 * Generates ",x" unless x==__ in which case it generates an empty string
 */
#define OPTIONS_ADD_CHAR(x) OPTIONS_PROCESS_CHAR((OPTIONS_NO_CHARACTER ## x, "," #x))

#define OPTIONS_NO_CHARACTER__ ~,) OPTIONS_DROP( ~
#define OPTIONS_DROP(x,y)
#define OPTIONS_PROCESS_CHAR(x) OPTIONS_PROCESS_CHAR2 x
#define OPTIONS_PROCESS_CHAR2(x,y) y

#define OPTIONS_OPTIONS_NONE(type) EMPTY
#define OPTIONS_OPTIONS_XEMPTY
#define OPTIONS_OPTIONS_D(def) D( def OPTIONS_OPTIONS_CLOSE
#define OPTIONS_OPTIONS_I(imp) I( imp OPTIONS_OPTIONS_CLOSE
#define OPTIONS_OPTIONS_DI(def,imp) DI( def, imp OPTIONS_OPTIONS_CLOSE
#define OPTIONS_OPTIONS_MULTI(type) MULTI
#define OPTIONS_OPTIONS_XMULTI ->multitoken()
#define OPTIONS_OPTIONS_CLOSE(type) ,type)
#define OPTIONS_OPTIONS_X(x) OPTIONS_OPTIONS_XX(OPTIONS_OPTIONS_ ## x)
#define OPTIONS_OPTIONS_XX(x) OPTIONS_OPTIONS_XXX(x)
#define OPTIONS_OPTIONS_XXX(x) OPTIONS_OPTIONS_X ## x
#define OPTIONS_OPTIONS_XD(def,type) ->default_value (SLOM::convert<type >(def),SLOM::convert_str(def))
#define OPTIONS_OPTIONS_XI(imp,type) ->implicit_value(SLOM::convert<type >(imp),SLOM::convert_str(imp))
#define OPTIONS_OPTIONS_XDI(def,imp,type) OPTIONS_OPTIONS_XD(def,type) OPTIONS_OPTIONS_XI(imp,type)

#define REGISTER_OPTION(type, name, letter, description, options) \
	(#name OPTIONS_ADD_CHAR(letter), SLOM::po::value(&name) OPTIONS_OPTIONS_X(options  (type)), description)


#define DECLARE_OPTION(type_, name, letter, description, options) \
		MTK::internal::UnalignedType<type_ >::type name; \
		self_& set_ ## name(const type_& value_) { name = value_; return *this; }


//REGISTER_OPTION(char, foo, f, "description", NONE )


#endif /* MAKEOPTIONS_HPP_ */
