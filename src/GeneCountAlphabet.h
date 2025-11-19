#ifndef GENE_COUNT_ALPHABET_H
#define GENE_COUNT_ALPHABET_H


#include <Bpp/Seq/Alphabet/AbstractAlphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetState.h>
#include <Bpp/Seq/Alphabet/AlphabetExceptions.h>
#include <Bpp/Text/TextTools.h>


namespace bpp
{

class GeneCountAlphabet :
  public AbstractAlphabet
{
private:
  unsigned int MIN_;
  unsigned int MAX_;
  

public:
  // class constructor
  GeneCountAlphabet(unsigned int max, unsigned int min = 0);

  GeneCountAlphabet(const GeneCountAlphabet& bia) : AbstractAlphabet(bia), MIN_(bia.MIN_), MAX_(bia.MAX_) {}

  GeneCountAlphabet& operator=(const GeneCountAlphabet& bia)
  {
    AbstractAlphabet::operator=(bia);
    MIN_=bia.MIN_;
    MAX_=bia.MAX_;

    return *this;
  }

  GeneCountAlphabet* clone() const
  {
    return new GeneCountAlphabet(*this);
  }
  // class destructor
  virtual ~GeneCountAlphabet() {}

public:
  unsigned int getSize() const { return MAX_ - MIN_ +1; }

  unsigned int getNumberOfTypes() const { return MAX_ - MIN_ + 1; }
  
  std::string getAlphabetType() const { return "Integer(MIN=" + TextTools::toString(MIN_) + ", MAX=" +  TextTools::toString(MAX_) +")"; }
  
  int getUnknownCharacterCode() const { return static_cast<int>(MAX_+1); }
  
  bool isUnresolved(int state) const { return state == static_cast<int>(MAX_+1); }
  
  bool isUnresolved(const std::string& state) const { return state == "X"; }
  
  unsigned int getMin() const { return MIN_; }

  unsigned int getMax() const { return MAX_; }

  bool isResolvedIn(int state1, int state2) const;
    
  std::vector<int> getAlias(int state) const;

  std::vector<std::string> getAlias(const std::string& state) const;
  


};
}
#endif // GENE_COUNT_ALPHABET_H