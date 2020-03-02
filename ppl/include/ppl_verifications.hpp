

#ifndef PPL_VERIFICATIONS_HPP
#define PPL_VERIFICATIONS_HPP


#include "ppl_headers.hpp"



namespace ppl
{

CONST ppl::LD eps{1e-17};

uint64_t prec_call(CONST ppl::LD& val)
{
	if(val>=1)
		return 1;
	return (uint64_t)std::ceil(-std::log10(val)-eps);
}


template<class T> inline CONST char* __Tn(void)
{
	if(strlen(typeid(T).name()) == 1)
	{
		switch(*typeid(T).name())
		{
			case 0x64: return "double";
			case 0x65: return "long double";
			case 0x66: return "float";
			default: return typeid(T).name();
		}

	}
	return typeid(T).name();
}






} // namespace ppl



#endif // PPL_VERIFICATIONS_H




