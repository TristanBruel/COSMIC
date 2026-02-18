***
      SUBROUTINE assign_remnant(zpars,mc,mcbagb,mass,kidx,caseMT,
     &                          mt,kw,bhspin)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      common /fall/fallback
      REAL*8 fallback
      REAL ran3,xx
      EXTERNAL ran3
      real*8 zpars(20)

      real*8 avar,bvar
      real*8 mc,mcbagb,mass,mt,mc_tot,met
      real*8 frac,kappa,sappa,alphap,polyfit
      real*8 m_proto,m_FeNi,m_fb,bhspin,mrem,mch,dMppi
      integer kw,kidx,caseMT
      real*8 logz
      real*8 Mco1,Mco2,Mco3,McoNS1,McoNS2


* Inputs
*       zpars      : Array of metallicity dependent parameters
*       mc         : CO core mass before SN
*       mcbagb     : Core mass at the base of the AGB
*       mass       : Previous epoch mass of the star
*       kidx       : Index of the star in the pisn track arrays
*       caseMT     : Mass transfer case of the star

* Outputs
*       mt         : Remnant mass after SN
*       kw         : Stellar (remnant) type
*       bhspin     : Dimensionless spin parameter of BH remnant

* total core mass before SN (CO + He layers)
      mc_tot = mc_co(kidx) + mc_he(kidx)

* Set the Chandrasekhar mass
      mch = 1.44d0 !set here owing to AIC ECSN model.

* Check if remnant is below Chandrasekhar mass
      if(mc_co(kidx).lt.mch)then
         if(ifflag.ge.1)then
*
* Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800.
*
            if(zpars(14).ge.1.0d-08)then
               mc = MIN(0.36d0+0.104d0*mass,0.58d0+0.061d0*mass)
               mc = MAX(0.54d0+0.042d0*mass,mc)
               if(mass.lt.1.d0) mc = 0.46d0
            else
               mc = MIN(0.29d0+0.178d0*mass,0.65d0+0.062d0*mass)
               mc = MAX(0.54d0+0.073d0*mass,mc)
            endif
            mc = MIN(mch,mc)
         endif

         mt = mc
         if(ecsn.gt.0.d0.and.mcbagb.lt.ecsn_mlow)then
            kw = 11
         elseif(ecsn.eq.0.d0.and.mcbagb.lt.1.6d0)then !double check what this should be. should be ecsn_mlow. Remember need to add option if ecsn = 0 (i.e. no ECSN!!!)
*
* Zero-age Carbon/Oxygen White Dwarf
*
            kw = 11
         elseif(ecsn.gt.0.d0.and.mcbagb.ge.ecsn_mlow.and.
     &          mcbagb.le.ecsn.and.mc.lt.1.08d0)then
            kw = 11
*               elseif(mcbagb.ge.1.6d0.and.mcbagb.le.2.5d0.and.
*                      mc.lt.1.08d0)then !can introduce this into code at some point.
*                  kw = 11

         else
*
* Zero-age Oxygen/Neon White Dwarf
*
            kw = 12
         endif
         mass = mt
*
      else
         if(ecsn.gt.0.d0.and.mcbagb.lt.ecsn_mlow)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
            kw = 15
         elseif(ecsn.eq.0.d0.and.mcbagb.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
            kw = 15
         else
*
* Use remnant mass given by Hurley+2000
*
            if(remnantflag.eq.0)then
               mt = 1.17d0 + 0.09d0*mc_co(kidx)
            elseif(remnantflag.eq.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407. Equation 1
*
*   TW: Belczynski+02 states that FeNi core mass comes from Woosley 1986 (Nucleosynthesis and Stellar Evolution)
*   As far I can tell, this is a linear fit to Table 6. But that table is He core mass, not CO core mass.
*   I also get a different slope, so I imagine there was a conversion to CO core mass somewhere in their fit?
               if(mc_co(kidx).lt.2.5d0)then
                  m_FeNi = 0.161767d0*mc_co(kidx) + 1.067055d0
               else
                  m_FeNi = 0.314154d0*mc_co(kidx) + 0.686088d0
               endif
               if(mc_co(kidx).le.5.d0)then
                  mt = m_FeNi
                  fallback = 0.d0
               elseif(mc_co(kidx).lt.7.6d0)then
*   TW: the fallback fraction is assumed to linearly increase from 0 to 1 between MCO=5 and MCO=7.6, hence the 2.6
                  fallback = (mc_co(kidx) - 5.d0)/2.6d0
                  mt = m_FeNi + fallback * (mt - m_FeNi)
               elseif(mc_co(kidx).gt.7.60)then
                  fallback = 1.d0
               endif
            elseif(remnantflag.eq.2)then
*
* Use NS/BH masses given by Belczynski+08. PK.
*
               ! calculate m_FeNi following Eq. 1 (fit to Timmes+1996)
               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                             mcbagb.ge.ecsn_mlow)then
                  m_FeNi = 1.38d0
               elseif(mc_co(kidx).lt.4.82d0)then
                  m_FeNi = 1.5d0
               elseif(mc_co(kidx).ge.4.82d0
     &                .and.mc_co(kidx).lt.6.31d0)then
                  m_FeNi = 2.11d0
               elseif(mc_co(kidx).ge.6.31d0
     &                .and.mc_co(kidx).lt.6.75d0)then
                  m_FeNi = 0.69*mc_co(kidx) - 2.26d0
               elseif(mc_co(kidx).ge.6.75d0)then
                  m_FeNi = 0.37*mc_co(kidx) - 0.07d0
               endif
               ! now calculate the remnant mass after fallback (Eq. 2)
               if(mc_co(kidx).le.5.d0)then
                  mt = m_FeNi
                  fallback = 0.d0
               elseif(mc_co(kidx).lt.7.6d0)then
                  fallback = (mc_co(kidx) - 5.d0) / 2.6d0
                  mt = m_FeNi + fallback * (mt - m_FeNi)
               elseif(mc_co(kidx).gt.7.60)then
                  fallback = 1.d0
               endif

            elseif(remnantflag.eq.3)then
*
* Use the "Rapid" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*              We use the updated proto-core mass from Giacobbo & Mapelli 2020
               m_proto = 1.1d0

               ! Calculate remnant mass from Eq. 16 + 17
               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                             mcbagb.ge.ecsn_mlow)then
                  mt = 1.38d0   ! ECSN fixed mass, no fallback
               elseif(mc_co(kidx).lt.2.5d0)then
                  fallback = 0.2d0 / (mt - m_proto)
                  mt = m_proto + 0.2d0
               elseif(mc_co(kidx).lt.6.d0)then
                  m_fb = 0.286d0 * mc_co(kidx) - 0.514d0
                  fallback = m_fb / (mt - m_proto)
                  mt = m_proto + m_fb
               elseif(mc_co(kidx).lt.7.d0)then
                  fallback = 1.d0
               elseif(mc_co(kidx).lt.11.d0)then
                  avar = 0.25d0 - (1.275 / (mt - m_proto))
                  bvar = 1.d0 - 11.d0*avar
                  fallback = avar*mc_co(kidx) + bvar
                  mt = m_proto + fallback*(mt - m_proto)
               elseif(mc_co(kidx).ge.11.d0)then
                  fallback = 1.d0
               endif
*              if the user requests it, limit the final remnant mass to
*              is the total **core** mass, not the total stellar mass
               if(fryer_mass_limit.eq.1)then
                  mt = min(mt, mc_tot)
               endif
            elseif(remnantflag.eq.4)then
*
* Use the "Delayed" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*              Calculate the proto-core mass following Eq. 18
               if(mc_co(kidx).lt.3.5d0)then
                  m_proto = 1.2d0
               elseif(mc_co(kidx).lt.6.d0)then
                  m_proto = 1.3d0
               elseif(mc_co(kidx).lt.11.d0)then
                  m_proto = 1.4d0
               elseif(mc_co(kidx).ge.11.d0)then
                  m_proto = 1.6d0
               endif

               ! Calculate remnant mass from Eq. 19 + 20
               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                             mcbagb.ge.ecsn_mlow)then
                  mt = 1.38d0   ! ECSN fixed mass, no fallback
               elseif(mc_co(kidx).lt.2.5d0)then
                  fallback = 0.2d0 / (mt - m_proto)
                  mt = m_proto + 0.2d0
               elseif(mc_co(kidx).lt.3.5d0)then
                  m_fb = 0.5d0 * mc_co(kidx) - 1.05d0
                  fallback = m_fb / (mt - m_proto)
                  mt = m_proto + m_fb
               elseif(mc_co(kidx).lt.11.d0)then
                  avar = 0.133d0 - (0.093d0 / (mt - m_proto))
                  bvar = 1.d0 - 11.d0*avar
                  fallback = avar*mc_co(kidx) + bvar
                  mt = m_proto + fallback*(mt - m_proto)
               elseif(mc_co(kidx).ge.11.d0)then
                  fallback = 1.d0
               endif
*              if the user requests it, limit the final remnant mass to
*              is the total **core** mass, not the total stellar mass
               if(fryer_mass_limit.eq.1)then
                  mt = min(mt, mc_tot)
               endif
            elseif(remnantflag.eq.5)then
               met = 10**(LOG10(zpars(14))/0.4)
*
* Use the Explodability Criteria from (Maltsev et al. 2025, A&A, 700,A20)
* with linear interpolation of the fallback fraction between direct BHs and NSs
*
*              Get the proto-core mass
               if(mc_co(kidx).le.3.5d0)then
                  m_proto = 1.2d0
               elseif(mc_co(kidx).le.6.d0)then
                  m_proto = 1.3d0
               elseif(mc_co(kidx).le.11.d0)then
                  m_proto = 1.4d0
               elseif(mc_co(kidx).gt.11.d0)then
                  m_proto = 1.6d0
               endif

*              Always Neutron Stars
               if(mc_co(kidx).lt.5.62d0)then
                  mt = mxns
*              Always Black Holes
               elseif(mc_co(kidx).gt.16.18d0)then
                  fallback = 1.d0
*              Identify the case of mass transfer
*              and compute the different ranges of Mco
               else
*                 Value of solar metallicity from Asplund et al. 2009
*                 extrapolate only between 1/20 and 1 [Zsun]
                  logz=MAX(log10(met/0.01432d0),log10(1.d0/20.d0))
*                 logz=MIN(logz,0.d0)
                  if(caseMT.eq.1)then
                     Mco1 = 7.4d0 + (7.4d0-6.9d0)*logz
                     Mco2 = 8.4d0 + (8.4d0-7.4d0)*logz
                     Mco3 = 15.4d0 + (15.4d0-13.7d0)*logz
                     McoNS1 = 11.1d0 + (11.1d0-10.4d0)*logz
                     McoNS2 = 12.1d0 + (12.1d0-11.1d0)*logz
                  elseif(caseMT.eq.2)then
                     Mco1 = 7.7d0 + (7.7d0-6.9d0)*logz
                     Mco2 = 8.3d0 + (8.3d0-7.9d0)*logz
                     Mco3 = 15.2d0 + (15.2d0-13.7d0)*logz
                     McoNS1 = 9.9d0 + (9.9d0-9.3d0)*logz
                     McoNS2 = 10.3d0 + (10.3d0-10.3d0)*logz
                  elseif(caseMT.eq.3)then
                     Mco1 = 6.6d0 + (6.6d0-6.3d0)*logz
                     Mco2 = 7.1d0 + (7.1d0-7.1d0)*logz
                     Mco3 = 13.2d0 + (13.2d0-12.3d0)*logz
                     McoNS1 = 9.6d0 + (9.6d0-8.9d0)*logz
                     McoNS2 = 11.7d0 + (11.7d0-9.5d0)*logz
                  else
                     Mco1 = 6.6d0 + (6.6d0-6.1d0)*logz
                     Mco2 = 7.2d0 + (7.2d0-6.6d0)*logz
                     Mco3 = 13.0d0 + (13.0d0-12.9d0)*logz
                     McoNS1 = 9.0d0 + (9.0d0-7.4d0)*logz
                     McoNS2 = 10.2d0 + (10.2d0-11.0d0)*logz
                  endif
*                 Range in which Mco lies:
*                 Neutron Stars
                  if(mc_co(kidx).lt.Mco1)then
                     mt = mxns
*                 Black Holes - direct collapse
                  elseif(mc_co(kidx).ge.Mco1.and. 
     &                             mc_co(kidx).lt.Mco2)then
                     fallback = 1.d0
*                 Neutron Stars
                  elseif(mc_co(kidx).ge.McoNS1.and.
     &                             mc_co(kidx).lt.McoNS2)then
                     mt = mxns
*                 Black Holes - direct collapse
                  elseif(mc_co(kidx).ge.Mco3)then
                     fallback = 1.d0
*                 Either Neutron Stars or Black Holes
                  elseif(mc_co(kidx).ge.Mco2.and.
     &                             mc_co(kidx).lt.McoNS1)then
                     xx = ran3(idum1)
*                    Neutron Stars
                     if(xx.gt.0.15d0)then
                        mt = mxns
*                    Black Holes
                     else
                        fallback = (mc_co(kidx)-McoNS1)/(Mco2-McoNS1)
                        fallback = 0.2d0+0.8d0*fallback
                        mt = m_proto + fallback*(mt - m_proto)
                     endif
                  else
                     xx = ran3(idum1)
*                    Neutron Stars
                     if(xx.gt.0.15d0)then
                        mt = mxns
*                    Black Holes
                     else
                        fallback = (mc_co(kidx)-McoNS2)/(Mco3-McoNS2)
                        fallback = 0.2d0+0.8d0*fallback
                        mt = m_proto + fallback*(mt - m_proto)
                     endif
                  endif
               endif
*              if the user requests it, limit the final remnant mass to
*              is the total **core** mass, not the total stellar mass
               if(fryer_mass_limit.eq.1)then
                  mt = min(mt, mc_tot)
               endif
            elseif(remnantflag.eq.6)then
               met = 10**(LOG10(zpars(14))/0.4)
*
* Model B from (Maltsev et al. 2025, A&A, 700,A20)
* with the Remnant Mass Relation from (Ugolini et al. 2025, A&A, 695,A122)
*

*              Always Neutron Stars
               if(mc_co(kidx).lt.5.62d0)then
                  mt = mxns
*              Always Black Holes
               elseif(mc_co(kidx).gt.16.18d0)then
                  fallback = 1.d0
*              Identify the case of mass transfer
*              and compute the different ranges of Mco
               else
*                 Value of solar metallicity from Asplund et al. 2009
*                 extrapolate only between 1/20 and 1 [Zsun]
                  logz=MAX(log10(met/0.01432d0),log10(1.d0/20.d0))
*                 logz=MIN(logz,0.d0)
                  if(caseMT.eq.1)then
                     Mco1 = 7.4d0 + (7.4d0-6.9d0)*logz
                     Mco2 = 8.4d0 + (8.4d0-7.4d0)*logz
                     Mco3 = 15.4d0 + (15.4d0-13.7d0)*logz
                  elseif(caseMT.eq.2)then
                     Mco1 = 7.7d0 + (7.7d0-6.9d0)*logz
                     Mco2 = 8.3d0 + (8.3d0-7.9d0)*logz
                     Mco3 = 15.2d0 + (15.2d0-13.7d0)*logz
                  elseif(caseMT.eq.3)then
                     Mco1 = 6.6d0 + (6.6d0-6.3d0)*logz
                     Mco2 = 7.1d0 + (7.1d0-7.1d0)*logz
                     Mco3 = 13.2d0 + (13.2d0-12.3d0)*logz
                  else
                     Mco1 = 6.6d0 + (6.6d0-6.1d0)*logz
                     Mco2 = 7.2d0 + (7.2d0-6.6d0)*logz
                     Mco3 = 13.0d0 + (13.0d0-12.9d0)*logz
                  endif
*                 Range in which Mco lies:
*                 Neutron Stars
                  if(mc_co(kidx).lt.Mco1)then
                     mt = mxns
*                 Black Holes - direct collapse
                  elseif(mc_co(kidx).ge.Mco1.and.
     &                             mc_co(kidx).le.Mco2)then
                     fallback = 1.d0
*                 Black Holes - direct collapse
                  elseif(mc_co(kidx).gt.Mco3)then
                     fallback = 1.d0
*                 Either Neutron Stars or Black Holes
                  else
                     xx = ran3(idum1)
                     fallback = 0.06d0*mc_co(kidx)-0.03d0
                     fallback=MAX(fallback,0.1d0)
*                    Neutron Stars
                     if(xx.gt.0.1d0)then
                        mt = mxns
*                    Black Holes
                     else
                        mt = MAX(fallback*mt,mxns+1.0d0)
                     endif
                  endif
               endif
*              if the user requests it, limit the final remnant mass to
*              is the total **core** mass, not the total stellar mass
               if(fryer_mass_limit.eq.1)then
                  mt = min(mt, mc_tot)
               endif
            endif

* Assign the BH spin based on the chosen prescription
            call assign_remnant_spin(mc, bhspin)

            ! convert from baryonic to gravitational mass
            call baryonic_to_gravitational_mass(mt, mrem)

* Determine whether a zero-age NS or BH is formed
            if(mrem.le.mxns)then
               mt = mrem
               mc = mt
               kw = 13
            else
               mt = mrem
               mc = mt
               kw = 14

* CLR - (Pulsational) Pair-Instability Supernova

* Belczynski+2016 prescription: just shrink any BH with a He core mass
* between 45 and 65 solar masses (provided the pisn flag is set at 45),
* and blow up anything between 65 and 135 solar masses.
* Cheap, but effective
               if(pisn.gt.0)then
                  if(mc_tot.ge.pisn.and.mc_tot.lt.65.d0)then
                     mt = pisn

                     ! convert from baryonic to gravitational mass
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

                     pisn_track(kidx)=6
                  elseif(mc_tot.ge.65.d0.and.mc_tot.lt.135.d0)then
                     mt = 0.d0
                     mc = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif
* The Spera+Mapelli2017 prescription is a tad more sophisticated:
* complex fitting formula to Stan Woosley's PSN models.  HOWEVER, these
* were done using the ZAMS mass/core mass/remnant mass relationships for
* SEVN, not BSE.  In other words, I woud be careful using this (and in
* practice, it doesn't vary that much from Belczynski's prescription,
* since the He core masses are the same in both)
               elseif(pisn.eq.-1)then
                  frac = mc_tot/mass
                  kappa = 0.67d0*frac + 0.1d0
                  sappa = 0.5226d0*frac - 0.52974d0
                  if(mc_tot.le.32.d0)then
                     alphap = 1.0d0
                  elseif(frac.lt.0.9d0.and.mc_tot.le.37.d0)then
                     alphap = 0.2d0*(kappa-1.d0)*mc_tot +
     &                        0.2d0*(37.d0 - 32.d0*kappa)
                     pisn_track(kidx)=6
                  elseif(frac.lt.0.9d0.and.mc_tot.le.60.d0)then
                     alphap = kappa
                     pisn_track(kidx)=6
                  elseif(frac.lt.0.9d0.and.mc_tot.lt.64.d0)then
                     alphap = kappa*(-0.25d0)*mc_tot+ kappa*16.d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.le.37.d0)then
                     alphap = sappa*(mc_tot - 32.d0) + 1.d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.le.56.d0.and.
     &                   sappa.lt.-0.034168d0)then
                     alphap = 5.d0*sappa + 1.d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.le.56.d0.and.
     &                   sappa.ge.-0.034168d0)then
                     alphap = (-0.1381d0*frac + 0.1309d0)*
     &                        (mc_tot - 56.d0) + 0.82916d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.lt.64.d0)then
                     alphap = -0.103645d0*mc_tot+ 6.63328d0
                     pisn_track(kidx)=6
                  elseif(mc_tot.ge.64.d0.and.mc_tot.lt.135.d0)then
                     alphap = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  elseif(mc_tot.ge.135.d0)then
                     alphap = 1.0d0
                  endif
                  mt = alphap*mt

                  ! convert from baryonic to gravitational mass
                  call baryonic_to_gravitational_mass(mt, mrem)
                  mt = mrem

* Fit (8th order polynomial) to Table 1 in Marchant+2018.
               elseif(pisn.eq.-2)then
                  if(mc_tot.ge.31.99d0.and.mc_tot.le.61.10d0)then
                     polyfit = -6.29429263d5
     &                      + 1.15957797d5*mc_tot
     &                      - 9.28332577d3*mc_tot**2d0
     &                      + 4.21856189d2*mc_tot**3d0
     &                      - 1.19019565d1*mc_tot**4d0
     &                      + 2.13499267d-1*mc_tot**5d0
     &                      - 2.37814255d-3*mc_tot**6d0
     &                      + 1.50408118d-5*mc_tot**7d0
     &                      - 4.13587235d-8*mc_tot**8d0
                     mt = polyfit
                     pisn_track(kidx)=6

                     ! convert from baryonic to gravitational mass
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

                  elseif(mc_tot.gt.61.10d0.and.
     &                   mc_tot.lt.124.12d0)then
                     mt = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif

* Fit (8th order polynomial) to Table 5 in Woosley2019.
               elseif(pisn.eq.-3)then
                  if(mc_tot.ge.29.53d0.and.mc_tot.le.60.12d0)then
                     polyfit = -3.14610870d5
     &                      + 6.13699616d4*mc_tot
     &                      - 5.19249710d3*mc_tot**2d0
     &                      + 2.48914888d2*mc_tot**3d0
     &                      - 7.39487537d0*mc_tot**4d0
     &                      + 1.39439936d-1*mc_tot**5d0
     &                      - 1.63012111d-3*mc_tot**6d0
     &                      + 1.08052344d-5*mc_tot**7d0
     &                      - 3.11019088d-8*mc_tot**8d0
                     mt = polyfit
                     pisn_track(kidx)=6

                     ! convert from baryonic to gravitational mass
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

                  elseif(mc_tot.gt.60.12d0.and.
     &                   mc_tot.lt.135.d0)then
                     mt = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif
* Apply the PPISN prescription from Renzo+2022 (https://ui.adsabs.harvard.edu/abs/2022RNAAS...6...25R/abstract)
* with the adaptations from Hendriks+2023 (https://scixplorer.org/abs/2023MNRAS.526.4130H/abstract)
* This is a top-down prescription, where we subtract mass from the total core mass
               elseif(pisn.eq.-4)then
                  if(mc_co(kidx).ge.38.d0+ppi_co_shift
     &               .and.mc_co(kidx).le.114.d0)then
*       Calculate DeltaM_PPI using Eq.6 from Hendriks+2023 (equivalently Eq.2 from Renzo+2022)
                     met = 10**(LOG10(zpars(14))/0.4)
                     dMppi = (0.0006d0 * LOG10(met) + 0.0054)
     &                      * (mc_co(kidx) - ppi_co_shift - 34.8d0)**3
     &                      - 0.0013 * (mc_co(kidx)
     &                                  - ppi_co_shift - 34.8d0)**2
     &                      + ppi_extra_ml
*       Set the remnant mass equal to the total core mass minus the PPI mass loss.
*       We use core mass not total mass because envelopes are expected to be removed by the first PPI pulse (e.g. Renzo+2020b)
                     mt = mc_tot - dMppi
                     
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

*       If the remnant mass is reduced below 10 Msun, assume a full PISN with no remnant
                     if(mt.lt.10.d0)then
                        mt = 0.0d0
                        kw = 15
                        pisn_track(kidx)=7
*       Otherwise we have a PPISN
                     else
                        pisn_track(kidx)=6
                     endif
*       For very large cores, we assume a full PISN with no remnant
                  elseif(mc_co(kidx).gt.114.d0)then
                     mt = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif
               endif

               mc = mt
* Store the initial BH mass for calculating the ISCO later
               if(Mbh_initial.eq.0)then
                  Mbh_initial = mt
               endif
            endif
         endif
      endif
*
      end


      SUBROUTINE baryonic_to_gravitational_mass(mt, mrem)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      real*8 mt, mrem

      ! remnantflag 0 and 1 already calculate gravitational mass
      if(remnantflag.le.1)then
         mrem = mt
      else
         ! negative values set the absolute maximum mass loss
         if(rembar_massloss.ge.0d0)then
            ! calculate Mrem from mt using quadratic formula
            ! mt - mrem = 0.075 mrem^2 (Lattimer & Yahil 1989, Timmes+1996)
            mrem = 6.6666667d0*(SQRT(1.d0 + 0.3d0*mt) - 1.d0)

            ! limit to maximum mass loss
            mrem = MAX(mrem, mt - rembar_massloss)

         ! positive values set the fractional mass loss
         else
            mrem = (1.d0 + rembar_massloss) * mt
         endif
      endif

      end


      SUBROUTINE assign_remnant_spin(mc, bhspin)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      real*8 ran3, mc, bhspin
      EXTERNAL ran3

* Set all BH spins equal to bhspinmag
      if(bhspinflag.eq.0)then
         bhspin = bhspinmag
* Randomly assign BH spins between 0 and bhspinmag
      elseif(bhspinflag.eq.1)then
         bhspin = ran3(idum1) * bhspinmag
* Assign BH spins based on Belczynski+17 prescription
      elseif(bhspinflag.eq.2)then
         if(mc.le.13.d0)then
            bhspin = 0.9d0
         elseif(mc.lt.27.d0)then
            bhspin = -0.064d0*mc + 1.736d0
         else
            bhspin = 0.0d0
         endif
      endif

      end
