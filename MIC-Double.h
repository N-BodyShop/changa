#ifndef __MIC_DOUBLE_H__
#define __MIC_DOUBLE_H__


#include <x86intrin.h>

#include<iostream>


class SSEDouble
{

   public: __m512d val;   

   public:
    
           SSEDouble() {} 
  
           SSEDouble(double d) { val = _mm512_set1_pd(d); }

           SSEDouble(double d0, double d1, double d2, double d3, double d4, double d5, double d6, double d7) {  val = _mm512_setr_pd(d0,d1,d2,d3,d4,d5,d6,d7); }

           /* Arithmetic Operators*/ 
           friend inline SSEDouble operator -(const SSEDouble &a) {SSEDouble c;c.val=_mm512_sub_pd(_mm512_setzero_pd(),a.val);return c;}

           friend inline SSEDouble operator +(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_add_pd(a.val,b.val);return c;}
                
           friend inline SSEDouble operator -(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_sub_pd(a.val,b.val);return c;}

           friend inline SSEDouble operator *(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_mul_pd(a.val,b.val);return c;}

           friend inline SSEDouble operator /(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_div_pd(a.val,b.val);return c;}

           friend inline SSEDouble sqrt      (const SSEDouble &a)                  { SSEDouble c;c.val= _mm512_sqrt_pd(a.val);return c;} 


          friend inline SSEDouble operator +(double a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_add_pd(_mm512_set1_pd(a),b.val);return c;}


          friend inline SSEDouble operator -(double a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_sub_pd(_mm512_set1_pd(a),b.val);return c;}

          friend inline SSEDouble operator *(double a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_mul_pd(_mm512_set1_pd(a),b.val);return c;}   
       
          friend inline SSEDouble operator /(double a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_div_pd(_mm512_set1_pd(a),b.val);return c;}

           inline SSEDouble& operator +=(const SSEDouble &a) {val= _mm512_add_pd(val,a.val);return *this;}
                
           inline SSEDouble& operator -=(const SSEDouble &a) {val= _mm512_sub_pd(val,a.val);return *this;}

           inline SSEDouble& operator *=(const SSEDouble &a) {val= _mm512_mul_pd(val,a.val);return *this;}

           inline SSEDouble& operator /=(const SSEDouble &a) {val= _mm512_div_pd(val,a.val);return *this;}

           /*Logical Operators*/

           friend inline SSEDouble operator &(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(a.val),_mm512_castpd_si512(b.val)));return c;}

           friend inline SSEDouble operator |(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_castsi512_pd(_mm512_or_epi64(_mm512_castpd_si512(a.val),_mm512_castpd_si512(b.val)));return c;}

           friend inline SSEDouble operator ^(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(a.val),_mm512_castpd_si512(b.val)));return c;}

           friend inline SSEDouble andnot (const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_castsi512_pd(_mm512_andnot_epi64(_mm512_castpd_si512(a.val),_mm512_castpd_si512(b.val)));return c;}

         /*Comparison Operators (not supported) */

            

            //friend inline SSEDouble operator <(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_cmplt_pd(a.val,b.val);return c;}
            friend inline SSEDouble operator <(const SSEDouble &a, const SSEDouble &b)
            {
              SSEDouble c;  // _mm512_castsi512_pd(_mm512_set1_epi64(0)),_mm512_castsi512_pd(_mm512_set1_epi64(-1))
              c.val= _mm512_mask_blend_pd(_mm512_cmp_pd_mask(a.val,b.val,_CMP_LT_OS),_mm512_set1_pd(0.0),_mm512_castsi512_pd(_mm512_set1_epi64(~0)));
              return c;
            }
            
            //friend inline SSEDouble operator >(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_cmpgt_pd(a.val,b.val);return c;}
            friend inline SSEDouble operator >(const SSEDouble &a, const SSEDouble &b)
            {
              SSEDouble c;
              c.val= _mm512_mask_blend_pd(_mm512_cmp_pd_mask(a.val,b.val,_CMP_GT_OS),_mm512_set1_pd(0.0),_mm512_castsi512_pd(_mm512_set1_epi64(~0)));
              return c;
            }
            //friend inline SSEDouble operator ==(const SSEDouble &a, const SSEDouble &b) {SSEDouble c;c.val= _mm512_cmpeq_pd(a.val,b.val);return c;}  
            friend inline SSEDouble operator ==(const SSEDouble &a, const SSEDouble &b)
            {
              SSEDouble c;
              c.val= _mm512_mask_blend_pd(_mm512_cmp_pd_mask(a.val,b.val,_CMP_EQ_OQ),_mm512_set1_pd(0.0),_mm512_castsi512_pd(_mm512_set1_epi64(~0)));
              return c;
            }
            //friend inline SSEDouble operator <(const SSEDouble &a, double b) {SSEDouble c;c.val= _mm512_cmplt_pd(a.val,_mm512_set1_pd(b));return c;} 
            friend inline SSEDouble operator <(const SSEDouble &a, const double b)
            {
              SSEDouble c;
              c.val= _mm512_mask_blend_pd(_mm512_cmp_pd_mask(a.val,_mm512_set1_pd(b),_CMP_LT_OS),_mm512_set1_pd(0.0),_mm512_castsi512_pd(_mm512_set1_epi64(~0)));
              return c;
            }
            //friend inline SSEDouble operator >(const SSEDouble &a, double b) {SSEDouble c;c.val= _mm512_cmpgt_pd(a.val,_mm512_set1_pd(b));return c;}
            friend inline SSEDouble operator >(const SSEDouble &a, const double b)
            {
              SSEDouble c;
              c.val= _mm512_mask_blend_pd(_mm512_cmp_pd_mask(a.val,_mm512_set1_pd(b),_CMP_GT_OS),_mm512_set1_pd(0.0),_mm512_castsi512_pd(_mm512_set1_epi64(~0)));
              return c;
            }
            
            friend inline SSEDouble max (const SSEDouble &a, SSEDouble &b) { SSEDouble c; c.val= _mm512_gmax_pd(a.val,b.val);return c;}
            

        /*Masking Operations */
                                                                                          /// Flip this order for opposite behavior
          friend inline int movemask( const SSEDouble &a) {return _mm512_mask2int(_mm512_cmpeq_epi64_mask(_mm512_castpd_si512(a.val),_mm512_set1_epi64(~0)));}


        /*Store Operations*/
            // should use _mm512_storeu_pd but doesn't exist 
        //  friend inline void storeu(double *p, const SSEDouble &a) { _mm512_storeu_pd(p,a.val);}
          friend inline void storeu(double *p, const SSEDouble &a) { _mm512_packstorelo_pd(p,a.val); _mm512_packstorehi_pd(p+16,a.val); }


          friend inline void store(double *p, const SSEDouble &a) { _mm512_store_pd(p,a.val);}

       //   void display();


 

};

#endif //__MIC_DOUBLE_H__
