#ifndef NUMERICAL_AND_MATHEMATICAL_CONSTANTS_H
#define NUMERICAL_AND_MATHEMATICAL_CONSTANTS_H

/// Provide numerical constants.
namespace constants {
template<class Number=double>
inline constexpr Number pi()
{
    return 3.14159265358979323846264338327950288419716939937510;
}

/// mass of uncharged pion
template<class Number=double>
inline constexpr Number mass_pi_0()
{return 0.134977;}

/// mass of charged pion
template<class Number=double>
inline constexpr Number mass_pi()
{ return 0.13957018;}

/// mass of eta meson
template<class Number=double>
inline constexpr Number mass_eta()
{return 0.547862;}

/// mass of eta' meson
template<class Number=double>
inline constexpr Number mass_etap()
{return 0.95778;}

/// mass of kaon
template<class Number=double>
inline constexpr Number mass_kaon()
{return 0.493677;}

/// mass of phi
template<class Number=double>
inline constexpr Number mass_phi()
{return 1.019461;}

/// mass of omega
template<class Number=double>
inline constexpr Number mass_omega()
{return 0.78266;}

} // constants


#endif // NUMERICAL_AND_MATHEMATICAL_CONSTANTS_H
