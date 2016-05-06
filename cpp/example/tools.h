#ifndef TOOLS_H_
#define TOOLS_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <stdio.h>



inline std::string make_filename( const std::string& basename, int index,
		const std::string& ext ) {
	std::ostringstream result;
	result << basename << std::setfill('0') << std::setw(3) << index << ext;
	return result.str();
}


struct nullstream : std::ostream {
	struct nullbuf: std::streambuf {
		int overflow(int c) { return traits_type::not_eof(c); }
	} m_sbuf;
	nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
};


/**
 * Very fast double to ASCII conversion, output will default to fixed precision,
 * but it will fall back for values too big.
 */
static inline
char* ftoa(double d, char * buf){
	if(d<0) {
		*buf++ = '-';
		d = -d;
	} else {
		*buf++ = ' ';
	}
	if(d >= (1e4 - 0.5e-5) ){
		int len = sprintf(buf,"%g", d);
		return buf+len;
	}
	
	d += 0.5e-5;
	d *= (1<<28)*1e-3;
	unsigned int z = d+.5;
	*buf++ = '0' + (z>>28);
	z &= ((1<<28)-1);
	for(int i=0; i<3; ++i){
		z*=10;
		*buf++ = '0' + (z>>28);
		z &= ((1<<28)-1);
	}
	*buf++ = '.';
	for(int i=0; i<5; ++i){
		z*=10;
		*buf++ = '0' + (z>>28);
		z &= ((1<<28)-1);
	}
	
	return buf;
}

static inline
char* itoa_half(unsigned int z, char *buf){
//	typedef unsigned int U;
	
	z = z * ((1<<28) / 10000 + 1) - (z >> 2);
	for(int i=0; i<5; ++i){
		*buf++ = '0' + (z>>28);
		z &= ((1<<28)-1);
		z*=10;
	}

	return buf;
}

static inline
char* itoa(unsigned int i, char *buf){
	
	if(i >= 100000){
		unsigned int u = i / 100000;
		i = i % 100000;
		buf = itoa_half(u, buf);
	}
	buf = itoa_half(i, buf);
	
	return buf;
}



#endif /* TOOLS_H_ */
