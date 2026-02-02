#include "GeneCountAlphabet.h"

using namespace bpp;

GeneCountAlphabet::GeneCountAlphabet(unsigned int max, unsigned int min) : MIN_(min), MAX_(max)
{
  registerState(new AlphabetState(-1, "-", "Gap"));

  for (int i = static_cast<int>(MIN_); i <= static_cast<int>(MAX_); ++i)
  {
    registerState(new AlphabetState(i, TextTools::toString(i), ""));
  }

  // M+ capped state
  registerState(new AlphabetState(MAX_ + 1, "M+", "Capped state"));

  // Optional unresolved state
  registerState(new AlphabetState(MAX_ + 2, "X", "Unresolved state"));
}
/********************************************************************************/
std::vector<int> GeneCountAlphabet::getAlias(int state) const{
  if (!isIntInAlphabet(state)) throw BadIntException(state, "GeneCountAlphabet::getAlias(int): Specified integer unknown.", this);
  std::vector<int> v;
  if (isUnresolved(state)){
    for (int i = static_cast<int>(MIN_); i <= static_cast<int>(MAX_+1); i++){
      v.push_back(i);
    }
  }else{
    v.push_back(state);
  }
  return v;
}

/********************************************************************************/
std::vector<std::string> GeneCountAlphabet::getAlias(const std::string& state) const{
  if (!isCharInAlphabet(state)) throw BadCharException(state, "GeneCountAlphabet::getAlias(char): Specified integer unknown.", this);
  std::vector<std::string> v;
  if (isUnresolved(state)){
    for (int i = static_cast<int>(MIN_); i <= static_cast<int>(MAX_+1); i++){
      v.push_back(intToChar(i));
    }
  }else if (charToInt(state) <= static_cast<int>(MAX_)){
    v.push_back(state);
  }else{
    throw Exception("GeneCountAlphabet::getAlias(char): unknown state!");

  } 
  return v;
}
/********************************************************************************/
bool GeneCountAlphabet::isResolvedIn(int state1, int state2) const{
  if (state1 < 0 || !isIntInAlphabet(state1))
    throw BadIntException(state1, "GeneCountAlphabet::isResolvedIn(int, int): Specified base unknown.", this);

  if (state2 < 0 || !isIntInAlphabet(state2))
    throw BadIntException(state2, "GeneCountAlphabet::isResolvedIn: Specified base unknown.", this);

  if (isUnresolved(state2))
    throw BadIntException(state2, "GeneCountAlphabet::isResolvedIn: Unresolved base.", this); 

  if (state2 > static_cast<int>(MAX_+1))
    throw IndexOutOfBoundsException("GeneCountAlphabet::isResolvedIn", state2, 0, getNumberOfTypes() - 1);
  std::vector<int> states = getAlias(state1);
  for (size_t j = 0; j < states.size(); j++)
  {
     if (state2 == states[j]){
        return true;

     }
  }
  return false;
}