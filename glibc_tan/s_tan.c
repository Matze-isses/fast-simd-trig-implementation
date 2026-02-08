/*
 * IBM Accurate Mathematical Library
 * written by International Business Machines Corp.
 * Copyright (C) 2001-2025 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, see <https://www.gnu.org/licenses/>.
 */
/*********************************************************************/
/*  MODULE_NAME: utan.c                                              */
/*                                                                   */
/*  FUNCTIONS: utan                                                  */
/*                                                                   */
/*  FILES NEEDED:dla.h endian.h mydefs.h utan.h                      */
/*               branred.c                                           */
/*               utan.tbl                                            */
/*                                                                   */
/*********************************************************************/

#include <stdint.h>
#include <errno.h>
#include <fenv.h>

/* If building standalone, these helpers/macros from glibc’s internal
   headers might be missing. Provide minimal replacements. */

/* HIGH_HALF / LOW_HALF: indices into mynumber.i[] for the 32-bit words
   of a double. The IBM code normally gets these from endian.h. */
#ifndef HIGH_HALF
# if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#  define HIGH_HALF 1
#  define LOW_HALF 0
# else
#  define HIGH_HALF 0
#  define LOW_HALF 1
# endif
#endif

/* glibc uses this to temporarily set rounding to nearest and restore it
   at function exit. For a simple standalone build, we ignore it. */
#ifndef SET_RESTORE_ROUND_53BIT
# define SET_RESTORE_ROUND_53BIT(mode) ((void)0)
#endif

/* glibc’s internal errno helper. */
#ifndef __set_errno
# define __set_errno(e) (errno = (e))
#endif

/* In glibc this is used to force underflow flags for tiny results.
   For testing the value of tan(x), it’s safe to make it a no-op. */
#ifndef math_check_force_underflow_nonneg
# define math_check_force_underflow_nonneg(x) ((void)0)
#endif
#include <endian.h>

#include "utan.tbl"
#include "dla.h"
#include "mydefs.h"
#include <math.h>
#include <stdio.h>
#include "bit_printing.h"










#ifndef SECTION
# define SECTION
#endif

/* tan with max ULP of ~0.619 based on random sampling.  */
double SECTION __tan (double x) {
#include "utan.h"
#include "utan.tbl"

  int ux, i, n;
  double a, da, a2, b, db, c, dc, fi, gi, pz,
	 s, sy, t, t1, t2, t3, t4, w, x2, xn, y, ya,
         yya, z0, z, z2;
  mynumber num, v;

  double retval;

  int __branred (double, double *, double *);

  SET_RESTORE_ROUND_53BIT (FE_TONEAREST);

  /* x=+-INF, x=NaN */
  num.d = x;
  ux = num.i[HIGH_HALF];

  if ((ux & 0x7ff00000) == 0x7ff00000)
    {
      if ((ux & 0x7fffffff) == 0x7ff00000)
	__set_errno (EDOM);
      retval = x - x;
      goto ret;
    }

  w = (x < 0.0) ? -x : x;

  /* (I) The case abs(x) <= 1.259e-8 */
  if (w <= g1.d)
    {
      math_check_force_underflow_nonneg (w);
      retval = x;
      goto ret;
    }

  /* (II) The case 1.259e-8 < abs(x) <= 0.0608 */
  if (w <= g2.d)
    {
      x2 = x * x;

      t2 = d9.d + x2 * d11.d;
      t2 = d7.d + x2 * t2;
      t2 = d5.d + x2 * t2;
      t2 = d3.d + x2 * t2;
      t2 *= x * x2;

      y = x + t2;
      retval = y;
      /* Max ULP is 0.504.  */
      goto ret;
    }

  





  // This is the interesting case
  
  /* (III) The case 0.0608 < abs(x) <= 0.787 */

  // w is exactly the input value.
  if (w <= g3.d)  // more precise g3.d = 0.78699970245361328
    {
      // mfftnhf = -15.5
      // mfftnhf = 1 10000000010 1111000000000000000000000000000000000000000000000000
      //
      // Due to the range max(i) = 185 and min(i) = 0. (int always rounding down)
      i = ((int) (mfftnhf.d + 256 * w)); 

      //To print out the values
      for (int j = 0; j < 201; j++) {
          printf("(%d, %#.17g),\n", j, xfg[j][0].d);
      }
      

      // xfg[i][0] does look like a linear increasing table, however it is not!                              
      // The values of z are in between -0.002 and 0.002 and jump wildly around it
      z = w - xfg[i][0].d;
      z2 = z * z;

      s = (x < 0.0) ? -1 : 1;

      // e0 is the 3th taylor coeff minus a small value 
      // e1 is the 5th taylor coeff plus a small value 
      // the first single z is for the 1-th coeff as it is 1
      // therefore pz is the taylor polynomial for the tan of degree 5 of z
      pz = z + z * z2 * (e0.d + z2 * e1.d);

      // Here it is hard to tell what the base function is. it looks really ugly.
      // Most likely there is no simple function which looks like this.
      fi = xfg[i][1].d;

      // gi are descrete developement points of the function: 
      //    Approx_gi(w) = 1 / w - 0.34851525202530603 * w 
      // This is shown in the approx_gi plot
      gi = xfg[i][2].d;

      print_double_bits(gi);

      t2 = pz * (gi + fi) / (gi - pz);

      // for (int j = 0; j < 201; j++) {
      //     printf("(%d, %#.17g),\n", j, xfg[j][2].d);
      // }

      // the values of fi almost behave like 1/x but as discrete number and faster decreasing.
      y = fi + t2;


      // simply adjusting the sign
      retval = (s * y);
      /* Max ULP is 0.60.  */
      goto ret;
    }







  /* (---) The case 0.787 < abs(x) <= 25 */
  if (w <= g4.d)
    {
      /* Range reduction by algorithm i */
      t = (x * hpinv.d + toint.d);
      xn = t - toint.d;
      v.d = t;
      t1 = (x - xn * mp1.d) - xn * mp2.d;
      n = v.i[LOW_HALF] & 0x00000001;
      da = xn * mp3.d;
      a = t1 - da;
      da = (t1 - a) - da;
      if (a < 0.0)
	{
	  ya = -a;
	  yya = -da;
	  sy = -1;
	}
      else
	{
	  ya = a;
	  yya = da;
	  sy = 1;
	}

      /* (VI) The case 0.787 < abs(x) <= 25,    0 < abs(y) <= 0.0608 */
  if (ya <= gy2.d) {
	  a2 = a * a;
	  t2 = d9.d + a2 * d11.d;
	  t2 = d7.d + a2 * t2;
	  t2 = d5.d + a2 * t2;
	  t2 = d3.d + a2 * t2;
	  t2 = da + a * a2 * t2;

	  if (n) {
	      /* -cot */
	      EADD (a, t2, b, db);
	      DIV2 (1.0, 0.0, b, db, c, dc, t1, t2, t3, t4);
	      y = c + dc;
	      retval = (-y);
	      /* Max ULP is 0.506.  */
	      goto ret;
	    } else {
	      /* tan */
	      y = a + t2;
	      retval = y;
	      /* Max ULP is 0.506.  */
	      goto ret;
	    }
  }

      /* (VII) The case 0.787 < abs(x) <= 25,    0.0608 < abs(y) <= 0.787 */

      i = ((int) (mfftnhf.d + 256 * ya));
      z = (z0 = (ya - xfg[i][0].d)) + yya;
      z2 = z * z;
      pz = z + z * z2 * (e0.d + z2 * e1.d);
      fi = xfg[i][1].d;
      gi = xfg[i][2].d;

      if (n) {
        /* -cot */
        t2 = pz * (fi + gi) / (fi + pz);
        y = gi - t2;
        retval = (-sy * y);
        /* Max ULP is 0.62.  */
        goto ret;
      } else {
        /* tan */
        t2 = pz * (gi + fi) / (gi - pz);
        y = fi + t2;
        retval = (sy * y);
        /* Max ULP is 0.62.  */
        goto ret;
      }
    }

  /* (---) The case 25 < abs(x) <= 1e8 */
  if (w <= g5.d)
    {
      /* Range reduction by algorithm ii */
      t = (x * hpinv.d + toint.d);
      xn = t - toint.d;
      v.d = t;
      t1 = (x - xn * mp1.d) - xn * mp2.d;
      n = v.i[LOW_HALF] & 0x00000001;
      da = xn * pp3.d;
      t = t1 - da;
      da = (t1 - t) - da;
      t1 = xn * pp4.d;
      a = t - t1;
      da = ((t - a) - t1) + da;
      EADD (a, da, t1, t2);
      a = t1;
      da = t2;
      if (a < 0.0)
	{
	  ya = -a;
	  yya = -da;
	  sy = -1;
	}
      else
	{
	  ya = a;
	  yya = da;
	  sy = 1;
	}

      /* (VIII) The case 25 < abs(x) <= 1e8,    0 < abs(y) <= 0.0608 */
      if (ya <= gy2.d)
	{
	  a2 = a * a;
	  t2 = d9.d + a2 * d11.d;
	  t2 = d7.d + a2 * t2;
	  t2 = d5.d + a2 * t2;
	  t2 = d3.d + a2 * t2;
	  t2 = da + a * a2 * t2;

	  if (n)
	    {
	      /* -cot */
	      EADD (a, t2, b, db);
	      DIV2 (1.0, 0.0, b, db, c, dc, t1, t2, t3, t4);
	      y = c + dc;
	      retval = (-y);
	      /* Max ULP is 0.506.  */
	      goto ret;
	    }
	  else
	    {
	      /* tan */
	      y = a + t2;
	      retval = y;
	      /* Max ULP is 0.506.  */
	      goto ret;
	    }
	}

      /* (IX) The case 25 < abs(x) <= 1e8,    0.0608 < abs(y) <= 0.787 */
      i = ((int) (mfftnhf.d + 256 * ya));
      z = (z0 = (ya - xfg[i][0].d)) + yya;
      z2 = z * z;
      pz = z + z * z2 * (e0.d + z2 * e1.d);
      fi = xfg[i][1].d;
      gi = xfg[i][2].d;

      if (n)
	{
	  /* -cot */
	  t2 = pz * (fi + gi) / (fi + pz);
	  y = gi - t2;
	  retval = (-sy * y);
	  /* Max ULP is 0.62.  */
	  goto ret;
	}
      else
	{
	  /* tan */
	  t2 = pz * (gi + fi) / (gi - pz);
	  y = fi + t2;
	  retval = (sy * y);
	  /* Max ULP is 0.62.  */
	  goto ret;
	}
    }

  /* (---) The case 1e8 < abs(x) < 2**1024 */
  /* Range reduction by algorithm iii */
  n = (__branred (x, &a, &da)) & 0x00000001;
  EADD (a, da, t1, t2);
  a = t1;
  da = t2;
  if (a < 0.0)
    {
      ya = -a;
      yya = -da;
      sy = -1;
    }
  else
    {
      ya = a;
      yya = da;
      sy = 1;
    }

  /* (X) The case 1e8 < abs(x) < 2**1024,    0 < abs(y) <= 0.0608 */
  if (ya <= gy2.d)
    {
      a2 = a * a;
      t2 = d9.d + a2 * d11.d;
      t2 = d7.d + a2 * t2;
      t2 = d5.d + a2 * t2;
      t2 = d3.d + a2 * t2;
      t2 = da + a * a2 * t2;
      if (n)
	{
	  /* -cot */
	  EADD (a, t2, b, db);
	  DIV2 (1.0, 0.0, b, db, c, dc, t1, t2, t3, t4);
	  y = c + dc;
	  retval = (-y);
	  /* Max ULP is 0.506.  */
	  goto ret;
	}
      else
	{
	  /* tan */
	  y = a + t2;
	  retval = y;
	  /* Max ULP is 0.507.  */
	  goto ret;
	}
    }

  /* (XI) The case 1e8 < abs(x) < 2**1024,    0.0608 < abs(y) <= 0.787 */
  i = ((int) (mfftnhf.d + 256 * ya));
  z = (z0 = (ya - xfg[i][0].d)) + yya;
  z2 = z * z;
  pz = z + z * z2 * (e0.d + z2 * e1.d);
  fi = xfg[i][1].d;
  gi = xfg[i][2].d;

  if (n)
    {
      /* -cot */
      t2 = pz * (fi + gi) / (fi + pz);
      y = gi - t2;
      retval = (-sy * y);
      /* Max ULP is 0.62.  */
      goto ret;
    }
  else
    {
      /* tan */
      t2 = pz * (gi + fi) / (gi - pz);
      y = fi + t2;
      retval = (sy * y);
      /* Max ULP is 0.62.  */
      goto ret;
    }

ret:
  return retval;
}

// gcc -O2 -Wall -Wextra main.c s_tan.c branred.c bit_printing.c -lm -o tan_test && ./tan_test 
