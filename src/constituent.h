#ifndef CONSTITUENT_H
#define CONSTITUENT_H

#include <array>
#include <string>
#include <unordered_map>

class Constituent {
 public:
  //  Constituent(const char *name) noexcept;
  //  Constituent(const std::string &name) noexcept;

  struct TC {
    std::string name;
    const double frequency;
    const std::array<double, 6> doodson;
    const double semi;
    const double isat;
    const int ishallow;
    const int nshallow;
    const double doodsonamp;
    const double earthreduc;
  };

  //  const char *name() const noexcept;
  //  bool isnull() const noexcept;

  //  double frequency() const;
  //  double doodson(size_t index) const;
  //  double semi() const;
  //  double isat() const;
  //  int ishallow() const;
  //  int nshallow() const;
  //  double doodsonamp() const;
  //  double earthreduc() const;

  template <typename T>
  static constexpr T nullvalue() {
    return -std::numeric_limits<T>::max();
  }

  static size_t constituentIndex(const std::string &name);
  static const std::array<std::array<int, 3>, 162> *deldood();
  static const std::array<double, 162> *phcorr();
  static const std::array<double, 162> *amprat();
  static const std::array<int, 162> *ilatfac();
  static const std::array<int, 162> *iconst();
  static std::vector<int> iconst_unique();
  static const std::array<const TC *, 146> *constituents();
  static const std::array<int, 251> *shallow_iconst();
  static const std::array<double, 251> *shallow_coef();
  static const std::array<int, 251> *shallow_iname();
  //  static const std::unordered_map<const char *, const TC *> *tideMap();

 private:
  //  static const TC *initialize(const std::string &name) noexcept;

  const TC *m_data;
};

#endif  // CONSTITUENT_H
