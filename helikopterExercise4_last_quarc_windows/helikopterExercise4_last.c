/*
 * helikopterExercise4_last.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterExercise4_last.mdl".
 *
 * Model version              : 1.71
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Thu Mar 19 12:38:17 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopterExercise4_last.h"
#include "helikopterExercise4_last_private.h"
#include <stdio.h>
#include "helikopterExercise4_last_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopterExercise4_las helikopterExercise4_last_B;

/* Continuous states */
ContinuousStates_helikopterExer helikopterExercise4_last_X;

/* Block states (auto storage) */
D_Work_helikopterExercise4_last helikopterExercise4_last_DWork;

/* Real-time model */
RT_MODEL_helikopterExercise4_la helikopterExercise4_last_M_;
RT_MODEL_helikopterExercise4_la *helikopterExercise4_last_M =
  &helikopterExercise4_last_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 5;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helikopterExercise4_last_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopterExercise4_last_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Gain2;
  real_T rtb_Sum[2];
  real_T rtb_Saturation_m;
  real_T rtb_Sum1_c[6];
  real_T rtb_Gain1;
  real_T rtb_degrad[6];
  real_T rtb_K_ed;
  int32_T tmp;
  real_T tmp_0[2];
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
    /* set solver stop time */
    if (!(helikopterExercise4_last_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopterExercise4_last_M->solverInfo,
                            ((helikopterExercise4_last_M->Timing.clockTickH0 + 1)
        * helikopterExercise4_last_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopterExercise4_last_M->solverInfo,
                            ((helikopterExercise4_last_M->Timing.clockTick0 + 1)
        * helikopterExercise4_last_M->Timing.stepSize0 +
        helikopterExercise4_last_M->Timing.clockTickH0 *
        helikopterExercise4_last_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopterExercise4_last_M)) {
    helikopterExercise4_last_M->Timing.t[0] = rtsiGetT
      (&helikopterExercise4_last_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise4_last_B.VandringLavpass =
    helikopterExercise4_last_P.VandringLavpass_C*
    helikopterExercise4_last_X.VandringLavpass_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
      if (!(++helikopterExercise4_last_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopterExercise4_last_DWork.ToFile_IWORK.Count*2)+1 < 100000000 )
      {
        FILE *fp = (FILE *) helikopterExercise4_last_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4_last_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopterExercise4_last_M->Timing.t[1];
          u[1] = helikopterExercise4_last_B.VandringLavpass;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4_last_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopterExercise4_last_DWork.ToFile_IWORK.Count)*2)+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file travel.mat.\n");
          }
        }
      }
    }

    /* S-Function (hil_read_encoder_block): '<S2>/HIL Read Encoder' */

    /* S-Function Block: helikopterExercise4_last/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder
        (helikopterExercise4_last_DWork.HILInitialize_Card,
         helikopterExercise4_last_P.HILReadEncoder_Channels, 3,
         &helikopterExercise4_last_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 =
          helikopterExercise4_last_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 =
          helikopterExercise4_last_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 =
          helikopterExercise4_last_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopterExercise4_last_B.KalibrerElev =
      helikopterExercise4_last_P.KalibrerElev_Gain * rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopterExercise4_last_B.Add = helikopterExercise4_last_B.KalibrerElev +
      helikopterExercise4_last_P.Constant_Value;

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterExercise4_last_B.KalibrerPitch =
      helikopterExercise4_last_P.KalibrerPitch_Gain * rtb_HILReadEncoder_o2;
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopterExercise4_last_X.Integrator_CSTATE >=
      helikopterExercise4_last_P.Integrator_UpperSat ) {
    helikopterExercise4_last_X.Integrator_CSTATE =
      helikopterExercise4_last_P.Integrator_UpperSat;
  } else if (helikopterExercise4_last_X.Integrator_CSTATE <=
             helikopterExercise4_last_P.Integrator_LowerSat ) {
    helikopterExercise4_last_X.Integrator_CSTATE =
      helikopterExercise4_last_P.Integrator_LowerSat;
  }

  rtb_Saturation_m = helikopterExercise4_last_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopterExercise4_last_B.KalibrerVandring =
      helikopterExercise4_last_P.KalibrerVandring_Gain * rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopterExercise4_last_P.VandringDeriv_D*
    helikopterExercise4_last_B.KalibrerVandring;
  rtb_VandringDeriv += helikopterExercise4_last_P.VandringDeriv_C*
    helikopterExercise4_last_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopterExercise4_last_P.TransferFcn4_D*
    helikopterExercise4_last_B.KalibrerPitch;
  rtb_Gain2 += helikopterExercise4_last_P.TransferFcn4_C*
    helikopterExercise4_last_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopterExercise4_last_P.TransferFcn5_D*
    helikopterExercise4_last_B.KalibrerElev;
  rtb_Gain1 += helikopterExercise4_last_P.TransferFcn5_C*
    helikopterExercise4_last_X.TransferFcn5_CSTATE;

  /* SignalConversion: '<Root>/TmpSignal ConversionAtdeg->radInport1' */
  rtb_Sum1_c[0] = helikopterExercise4_last_B.VandringLavpass;
  rtb_Sum1_c[1] = rtb_VandringDeriv;
  rtb_Sum1_c[2] = helikopterExercise4_last_B.KalibrerPitch;
  rtb_Sum1_c[3] = rtb_Gain2;
  rtb_Sum1_c[4] = helikopterExercise4_last_B.Add;
  rtb_Sum1_c[5] = rtb_Gain1;

  /* Gain: '<Root>/deg->rad' */
  for (tmp = 0; tmp < 6; tmp++) {
    rtb_degrad[tmp] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_degrad[tmp] += helikopterExercise4_last_P.degrad_Gain[6 * tmp_1 + tmp]
        * rtb_Sum1_c[tmp_1];
    }
  }

  /* Gain: '<S1>/K_ed' */
  rtb_K_ed = helikopterExercise4_last_P.K_ed_Gain * rtb_degrad[5];

  /* FromWorkspace: '<Root>/u_star' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise4_last_DWork.u_star_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise4_last_DWork.u_star_PWORK.TimePtr;
    int_T currTimeIndex = helikopterExercise4_last_DWork.u_star_IWORK.PrevIndex;
    real_T t = helikopterExercise4_last_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[80]) {
      currTimeIndex = 79;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helikopterExercise4_last_DWork.u_star_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Sum[0] = pDataValues[currTimeIndex];
          pDataValues += 81;
          rtb_Sum[1] = pDataValues[currTimeIndex];
          pDataValues += 81;
        } else {
          rtb_Sum[0] = pDataValues[currTimeIndex + 1];
          pDataValues += 81;
          rtb_Sum[1] = pDataValues[currTimeIndex + 1];
          pDataValues += 81;
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum[0] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 81;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum[1] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 81;
      }
    }
  }

  /* FromWorkspace: '<Root>/x_star' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise4_last_DWork.x_star_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise4_last_DWork.x_star_PWORK.TimePtr;
    int_T currTimeIndex = helikopterExercise4_last_DWork.x_star_IWORK.PrevIndex;
    real_T t = helikopterExercise4_last_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[80]) {
      currTimeIndex = 79;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helikopterExercise4_last_DWork.x_star_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T i1;
            real_T *y0 = rtb_Sum1_c;
            for (i1=0; i1 < 6; i1++) {
              y0[i1] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T i1;
            real_T *y0 = rtb_Sum1_c;
            for (i1=0; i1 < 6; i1++) {
              y0[i1] = pDataValues[currTimeIndex + 1];
              pDataValues += 81;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T i1;
          real_T *y0 = rtb_Sum1_c;
          for (i1=0; i1 < 6; i1++) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            y0[i1] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  /* Sum: '<S3>/Sum1' */
  for (tmp = 0; tmp < 6; tmp++) {
    rtb_Sum1_c[tmp] = rtb_degrad[tmp] - rtb_Sum1_c[tmp];
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/Gain'
   */
  for (tmp = 0; tmp < 2; tmp++) {
    tmp_0[tmp] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      tmp_0[tmp] += helikopterExercise4_last_P.Gain_Gain[(tmp_1 << 1) + tmp] *
        rtb_Sum1_c[tmp_1];
    }

    rtb_Sum[tmp] -= tmp_0[tmp];
  }

  /* Sum: '<S1>/Sum' */
  rtb_Gain1 = rtb_Sum[1] - rtb_degrad[4];

  /* Gain: '<S1>/K_ei' */
  helikopterExercise4_last_B.K_ei = helikopterExercise4_last_P.K_ei_Gain *
    rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_m = (helikopterExercise4_last_P.K_ep_Gain * rtb_Gain1 +
                      rtb_Saturation_m) + rtb_K_ed;

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_m = rt_SATURATE(rtb_Saturation_m,
    helikopterExercise4_last_P.Saturation_LowerSat,
    helikopterExercise4_last_P.Saturation_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
  }

  /* Sum: '<S4>/Sum' incorporates:
   *  Gain: '<S4>/K_pd'
   *  Gain: '<S4>/K_pp'
   *  Saturate: '<S4>/Saturation'
   *  Sum: '<S4>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Sum[0],
    helikopterExercise4_last_P.Saturation_LowerSat_c,
    helikopterExercise4_last_P.Saturation_UpperSat_d) - rtb_degrad[2]) *
    helikopterExercise4_last_P.K_pp_Gain - helikopterExercise4_last_P.K_pd_Gain *
    rtb_degrad[3];

  /* Gain: '<S5>/Gain2' incorporates:
   *  Sum: '<S5>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_m - rtb_Gain1) *
    helikopterExercise4_last_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopterExercise4_last_B.SatB = rt_SATURATE(rtb_Gain2,
    helikopterExercise4_last_P.SatB_LowerSat,
    helikopterExercise4_last_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
  }

  /* Gain: '<S5>/Gain1' incorporates:
   *  Sum: '<S5>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_m) *
    helikopterExercise4_last_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopterExercise4_last_B.Sat = rt_SATURATE(rtb_Gain1,
    helikopterExercise4_last_P.Sat_LowerSat,
    helikopterExercise4_last_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopterExercise4_last/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopterExercise4_last_DWork.HILWriteAnalog_Buffer[0] =
        helikopterExercise4_last_B.SatB;
      helikopterExercise4_last_DWork.HILWriteAnalog_Buffer[1] =
        helikopterExercise4_last_B.Sat;
      result = hil_write_analog
        (helikopterExercise4_last_DWork.HILInitialize_Card,
         helikopterExercise4_last_P.HILWriteAnalog_Channels, 2,
         &helikopterExercise4_last_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopterExercise4_last_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
    rt_ertODEUpdateContinuousStates(&helikopterExercise4_last_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helikopterExercise4_last_M->Timing.clockTick0)) {
    ++helikopterExercise4_last_M->Timing.clockTickH0;
  }

  helikopterExercise4_last_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopterExercise4_last_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopterExercise4_last_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopterExercise4_last_M->Timing.clockTick1)) {
      ++helikopterExercise4_last_M->Timing.clockTickH1;
    }

    helikopterExercise4_last_M->Timing.t[1] =
      helikopterExercise4_last_M->Timing.clockTick1 *
      helikopterExercise4_last_M->Timing.stepSize1 +
      helikopterExercise4_last_M->Timing.clockTickH1 *
      helikopterExercise4_last_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopterExercise4_last_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->VandringLavpass_CSTATE =
      helikopterExercise4_last_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->VandringLavpass_CSTATE +=
      (helikopterExercise4_last_P.VandringLavpass_A)*
      helikopterExercise4_last_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopterExercise4_last_X.Integrator_CSTATE <=
            helikopterExercise4_last_P.Integrator_LowerSat );
    usat = ( helikopterExercise4_last_X.Integrator_CSTATE >=
            helikopterExercise4_last_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopterExercise4_last_B.K_ei > 0)) ||
        (usat && (helikopterExercise4_last_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopterExer *)
        helikopterExercise4_last_M->ModelData.derivs)->Integrator_CSTATE =
        helikopterExercise4_last_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopterExer *)
        helikopterExercise4_last_M->ModelData.derivs)->Integrator_CSTATE = 0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->VandringDeriv_CSTATE =
      helikopterExercise4_last_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->VandringDeriv_CSTATE +=
      (helikopterExercise4_last_P.VandringDeriv_A)*
      helikopterExercise4_last_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->TransferFcn4_CSTATE =
      helikopterExercise4_last_B.KalibrerPitch;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->TransferFcn4_CSTATE +=
      (helikopterExercise4_last_P.TransferFcn4_A)*
      helikopterExercise4_last_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->TransferFcn5_CSTATE =
      helikopterExercise4_last_B.KalibrerElev;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_last_M->ModelData.derivs)->TransferFcn5_CSTATE +=
      (helikopterExercise4_last_P.TransferFcn5_A)*
      helikopterExercise4_last_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopterExercise4_last_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopterExercise4_last_P.Integrator_UpperSat = rtInf;
  helikopterExercise4_last_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopterExercise4_last_M, 0,
                sizeof(RT_MODEL_helikopterExercise4_la));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopterExercise4_last_M->solverInfo,
                          &helikopterExercise4_last_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopterExercise4_last_M->solverInfo, &rtmGetTPtr
                (helikopterExercise4_last_M));
    rtsiSetStepSizePtr(&helikopterExercise4_last_M->solverInfo,
                       &helikopterExercise4_last_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopterExercise4_last_M->solverInfo,
                 &helikopterExercise4_last_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopterExercise4_last_M->solverInfo,
                         &helikopterExercise4_last_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopterExercise4_last_M->solverInfo,
      &helikopterExercise4_last_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopterExercise4_last_M->solverInfo,
                          (&rtmGetErrorStatus(helikopterExercise4_last_M)));
    rtsiSetRTModelPtr(&helikopterExercise4_last_M->solverInfo,
                      helikopterExercise4_last_M);
  }

  rtsiSetSimTimeStep(&helikopterExercise4_last_M->solverInfo, MAJOR_TIME_STEP);
  helikopterExercise4_last_M->ModelData.intgData.f[0] =
    helikopterExercise4_last_M->ModelData.odeF[0];
  helikopterExercise4_last_M->ModelData.contStates = ((real_T *)
    &helikopterExercise4_last_X);
  rtsiSetSolverData(&helikopterExercise4_last_M->solverInfo, (void *)
                    &helikopterExercise4_last_M->ModelData.intgData);
  rtsiSetSolverName(&helikopterExercise4_last_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopterExercise4_last_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopterExercise4_last_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopterExercise4_last_M->Timing.sampleTimes =
      (&helikopterExercise4_last_M->Timing.sampleTimesArray[0]);
    helikopterExercise4_last_M->Timing.offsetTimes =
      (&helikopterExercise4_last_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopterExercise4_last_M->Timing.sampleTimes[0] = (0.0);
    helikopterExercise4_last_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopterExercise4_last_M->Timing.offsetTimes[0] = (0.0);
    helikopterExercise4_last_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopterExercise4_last_M,
             &helikopterExercise4_last_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopterExercise4_last_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopterExercise4_last_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopterExercise4_last_M, -1);
  helikopterExercise4_last_M->Timing.stepSize0 = 0.001;
  helikopterExercise4_last_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopterExercise4_last_M->Sizes.checksums[0] = (403469783U);
  helikopterExercise4_last_M->Sizes.checksums[1] = (471337330U);
  helikopterExercise4_last_M->Sizes.checksums[2] = (1547134978U);
  helikopterExercise4_last_M->Sizes.checksums[3] = (296910212U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopterExercise4_last_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopterExercise4_last_M->extModeInfo,
      &helikopterExercise4_last_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopterExercise4_last_M->extModeInfo,
                        helikopterExercise4_last_M->Sizes.checksums);
    rteiSetTPtr(helikopterExercise4_last_M->extModeInfo, rtmGetTPtr
                (helikopterExercise4_last_M));
  }

  helikopterExercise4_last_M->solverInfoPtr =
    (&helikopterExercise4_last_M->solverInfo);
  helikopterExercise4_last_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopterExercise4_last_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopterExercise4_last_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopterExercise4_last_M->ModelData.blockIO = ((void *)
    &helikopterExercise4_last_B);

  {
    helikopterExercise4_last_B.VandringLavpass = 0.0;
    helikopterExercise4_last_B.KalibrerElev = 0.0;
    helikopterExercise4_last_B.Add = 0.0;
    helikopterExercise4_last_B.KalibrerPitch = 0.0;
    helikopterExercise4_last_B.KalibrerVandring = 0.0;
    helikopterExercise4_last_B.K_ei = 0.0;
    helikopterExercise4_last_B.SatB = 0.0;
    helikopterExercise4_last_B.Sat = 0.0;
  }

  /* parameters */
  helikopterExercise4_last_M->ModelData.defaultParam = ((real_T *)
    &helikopterExercise4_last_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopterExercise4_last_X;
    helikopterExercise4_last_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopterExercise4_last_X, 0,
                  sizeof(ContinuousStates_helikopterExer));
  }

  /* states (dwork) */
  helikopterExercise4_last_M->Work.dwork = ((void *)
    &helikopterExercise4_last_DWork);
  (void) memset((void *)&helikopterExercise4_last_DWork, 0,
                sizeof(D_Work_helikopterExercise4_last));

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise4_last_DWork.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise4_last_DWork.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  helikopterExercise4_last_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopterExercise4_last_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopterExercise4_last_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopterExercise4_last_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopterExercise4_last_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 15;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }
}

/* Model terminate function */
void helikopterExercise4_last_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise4_last/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopterExercise4_last_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopterExercise4_last_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopterExercise4_last_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopterExercise4_last_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopterExercise4_last_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopterExercise4_last_P.HILInitialize_AOExit && is_switching)) {
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise4_last_P.HILInitialize_AOFinal;
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise4_last_P.HILInitialize_AOFinal;
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise4_last_P.HILInitialize_AOFinal;
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise4_last_P.HILInitialize_AOFinal;
      result = hil_write_analog
        (helikopterExercise4_last_DWork.HILInitialize_Card,
         &helikopterExercise4_last_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4_last_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
      }
    }

    hil_close(helikopterExercise4_last_DWork.HILInitialize_Card);
    helikopterExercise4_last_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopterExercise4_last_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_last_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4_last_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4_last_DWork.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4_last_M,
                          "Error writing header for ans to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_last_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      helikopterExercise4_last_DWork.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of GRT compatible call interface                                 *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helikopterExercise4_last_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopterExercise4_last_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopterExercise4_last_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopterExercise4_last_M->Sizes.numY = (0);/* Number of model outputs */
  helikopterExercise4_last_M->Sizes.numU = (0);/* Number of model inputs */
  helikopterExercise4_last_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopterExercise4_last_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopterExercise4_last_M->Sizes.numBlocks = (47);/* Number of blocks */
  helikopterExercise4_last_M->Sizes.numBlockIO = (8);/* Number of block outputs */
  helikopterExercise4_last_M->Sizes.numBlockPrms = (160);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise4_last_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopterExercise4_last_X.Integrator_CSTATE =
    helikopterExercise4_last_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopterExercise4_last_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopterExercise4_last_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopterExercise4_last_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise4_last/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("sensoray_model_626", "0",
                      &helikopterExercise4_last_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopterExercise4_last_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopterExercise4_last_P.HILInitialize_CKPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_ClockModes =
          &helikopterExercise4_last_DWork.HILInitialize_ClockModes[0];
        for (i1=0; i1 < 6; i1++) {
          dw_ClockModes[i1] = helikopterExercise4_last_P.HILInitialize_CKModes;
        }
      }

      result = hil_set_clock_mode
        (helikopterExercise4_last_DWork.HILInitialize_Card, (t_clock *)
         &helikopterExercise4_last_P.HILInitialize_CKChannels[0], 6U,
         (t_clock_mode *)
         &helikopterExercise4_last_DWork.HILInitialize_ClockModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear
      (helikopterExercise4_last_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
      return;
    }

    if ((helikopterExercise4_last_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopterExercise4_last_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helikopterExercise4_last_DWork.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMinimums[i1] = helikopterExercise4_last_P.HILInitialize_AILow;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helikopterExercise4_last_DWork.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMaximums[i1] = helikopterExercise4_last_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helikopterExercise4_last_DWork.HILInitialize_Card,
         &helikopterExercise4_last_P.HILInitialize_AIChannels[0], 16U,
         &helikopterExercise4_last_DWork.HILInitialize_AIMinimums[0],
         &helikopterExercise4_last_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise4_last_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopterExercise4_last_P.HILInitialize_AOPEnter && is_switching)) {
      helikopterExercise4_last_DWork.HILInitialize_AOMinimums[0] =
        helikopterExercise4_last_P.HILInitialize_AOLow;
      helikopterExercise4_last_DWork.HILInitialize_AOMinimums[1] =
        helikopterExercise4_last_P.HILInitialize_AOLow;
      helikopterExercise4_last_DWork.HILInitialize_AOMinimums[2] =
        helikopterExercise4_last_P.HILInitialize_AOLow;
      helikopterExercise4_last_DWork.HILInitialize_AOMinimums[3] =
        helikopterExercise4_last_P.HILInitialize_AOLow;
      helikopterExercise4_last_DWork.HILInitialize_AOMaximums[0] =
        helikopterExercise4_last_P.HILInitialize_AOHigh;
      helikopterExercise4_last_DWork.HILInitialize_AOMaximums[1] =
        helikopterExercise4_last_P.HILInitialize_AOHigh;
      helikopterExercise4_last_DWork.HILInitialize_AOMaximums[2] =
        helikopterExercise4_last_P.HILInitialize_AOHigh;
      helikopterExercise4_last_DWork.HILInitialize_AOMaximums[3] =
        helikopterExercise4_last_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopterExercise4_last_DWork.HILInitialize_Card,
         &helikopterExercise4_last_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4_last_DWork.HILInitialize_AOMinimums[0],
         &helikopterExercise4_last_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise4_last_P.HILInitialize_AOStart && !is_switching) ||
        (helikopterExercise4_last_P.HILInitialize_AOEnter && is_switching)) {
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise4_last_P.HILInitialize_AOInitial;
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise4_last_P.HILInitialize_AOInitial;
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise4_last_P.HILInitialize_AOInitial;
      helikopterExercise4_last_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise4_last_P.HILInitialize_AOInitial;
      result = hil_write_analog
        (helikopterExercise4_last_DWork.HILInitialize_Card,
         &helikopterExercise4_last_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4_last_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_last_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4_last_M,
                        "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4_last_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopterExercise4_last_DWork.ToFile_IWORK.Count = 0;
    helikopterExercise4_last_DWork.ToFile_IWORK.Decimation = -1;
    helikopterExercise4_last_DWork.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/u_star' */
  {
    static real_T pTimeValues[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.6455121993179131E-001, 8.2275609965895641E-002, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 4.6508511855277218E-001,
      4.0200755881036959E-001, 2.2808272411722139E-001, -1.7900893482856929E-002,
      -1.9359397850806631E-001, -2.8742275752476232E-001,
      -3.5036404535357996E-001, -4.0873729994601082E-001,
      -4.5421193253605296E-001, -4.6428084889936777E-001,
      -4.5496662634464108E-001, -4.5011378578512667E-001,
      -4.3314400202376824E-001, -3.9235017692082019E-001,
      -3.4742097717751153E-001, -3.1059903001704536E-001,
      -2.7192446675800230E-001, -2.2427589229108094E-001,
      -1.7511635158849584E-001, -1.3374516637370037E-001,
      -9.9631164949389464E-002, -6.7186794964896338E-002,
      -3.6375832819895169E-002, -1.1872048027525608E-002,
      1.2835456622208228E-003, 3.1512084577829680E-004, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5719423857623291E-003,
      1.2221744528993380E-002, 1.6638702459657106E-002, 2.3031297902052886E-002,
      3.1832678266465832E-002, 4.3668613460816942E-002, 5.9537225181344577E-002,
      8.0752100009053890E-002, 1.0825027178731837E-001, 1.4304083525633746E-001,
      1.8475954244149273E-001, 2.3060586368680336E-001, 2.7151577427895607E-001,
      2.8607575066473145E-001, 2.2748545611083540E-001, -2.6620032527911286E-005,
      -6.0600428787701833E-005, -4.0772964153704746E-007,
      2.7249278319789954E-005, -4.2782725177049062E-007, 2.7888082637255269E-005,
      -1.9446811047979072E-005, 4.3308101406108499E-005,
      -1.0781190692805759E-005, -1.4954024093787342E-005,
      -2.5941138801716148E-005, 8.6711963123805999E-005,
      -4.4567225445025789E-005, -8.1867139683634350E-005,
      1.1912315173645465E-004, -4.8517632751315844E-005,
      -2.3914462369836241E-005, 2.1749700778017312E-005, 1.1489590896809077E-005,
      -8.0525919108974698E-006, -8.5885390671500218E-006,
      2.0230202911757141E-005, -5.4992478255153272E-005, 2.2272791770139828E-005,
      -1.5122122644725121E-019, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopterExercise4_last_DWork.u_star_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise4_last_DWork.u_star_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise4_last_DWork.u_star_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/x_star' */
  {
    static real_T pTimeValues[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues[] = { 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1296829880716439E+000,
      3.1058636570353455E+000, 3.0701346604808979E+000, 3.0224959984083011E+000,
      2.9629476708175551E+000, 2.8914896777086598E+000, 2.8081220190816154E+000,
      2.7128446949364218E+000, 2.6056577052730789E+000, 2.4865610500915869E+000,
      2.3555547293919457E+000, 2.2126387431741552E+000, 2.0578130914382156E+000,
      1.8910777741841269E+000, 1.7166678136982314E+000, 1.5412660633833590E+000,
      1.3694437990094785E+000, 1.2033668452469342E+000, 1.0445741182378074E+000,
      8.9481738888200424E-001, 7.5540067887717921E-001, 6.2650264407886536E-001,
      5.0800655304943299E-001, 4.0010364294941236E-001, 3.0252847160467600E-001,
      2.1428419126287362E-001, 1.3435607274198039E-001, 6.1957757430275644E-002,
      -3.8342749023816927E-003, -6.4037834958831419E-002,
      -1.1948410003114680E-001, -1.7130755433462788E-001,
      -2.2140424138687309E-001, -2.7067382604442408E-001,
      -3.1695004599561255E-001, -3.5908510117458847E-001,
      -4.0454495898279330E-001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, -4.7638662072596753E-002, -9.5277324145193507E-002,
      -1.4291598621779025E-001, -1.9055464829038701E-001,
      -2.3819331036298377E-001, -2.8583197243558051E-001,
      -3.3347063450817727E-001, -3.8110929658077403E-001,
      -4.2874795865337079E-001, -4.7638662072596755E-001,
      -5.2402528279856431E-001, -5.7166394487116101E-001,
      -6.1930260694375772E-001, -6.6694126901635442E-001,
      -6.9763984194358275E-001, -7.0160700125949005E-001,
      -6.8728905749552127E-001, -6.6430781505017766E-001,
      -6.3517090803650667E-001, -5.9902691742321346E-001,
      -5.5766684001930023E-001, -5.1559213919325519E-001,
      -4.7398436411772937E-001, -4.3161164040008254E-001,
      -3.9030068537894552E-001, -3.5297712136720949E-001,
      -3.1971247408357284E-001, -2.8959326124681900E-001,
      -2.6316812933062933E-001, -2.4081424022579889E-001,
      -2.2178506028926154E-001, -2.0729381721392437E-001,
      -2.0038674820898078E-001, -1.9707833863020385E-001,
      -1.8510487980475407E-001, -1.6854022071590374E-001,
      -1.8183943123281943E-001, -2.2947809330541619E-001, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 3.3740945899817737E-001, 4.3603234642630842E-002,
      -1.5736919337648114E-001, -2.5258791667516917E-001,
      -3.2024511548688800E-001, -3.9726030091938530E-001,
      -4.5459055618182631E-001, -4.6244501582789826E-001,
      -4.5731301294136462E-001, -4.6572059944746313E-001,
      -4.5405064976218590E-001, -4.1022504762473611E-001,
      -3.6561330294874977E-001, -3.3104168499267028E-001,
      -2.9043986783792181E-001, -2.4569264660105017E-001,
      -2.0915061174948679E-001, -1.5927393425913389E-001,
      -7.5915920318042951E-002, -3.6362885296515123E-002,
      -1.3160084914073245E-001, -1.8206294718999652E-001,
      1.4617224834042220E-001, 5.2359877559829882E-001, -1.2614751684548914E-001,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      2.0943951023931953E+000, -1.5672700220995666E-016,
      -1.2881911772658114E-016, 1.3846294180621298E-016,
      -1.6301663779137805E-016, 8.0855223532245701E-017,
      -9.5637924689842119E-017, -8.2136767624558159E-017,
      -7.7942558326873029E-018, -5.0433882249827140E-017,
      3.0084062275783441E-016, -1.3901595072570175E-017,
      -1.4388527920375888E-016, -1.7197986584992306E-016,
      -7.4475726640048578E-001, -1.1752248974221862E+000,
      -8.0388971207644799E-001, -3.8087489319475204E-001,
      -2.7062879524687539E-001, -3.0806074172998910E-001,
      -2.2932102104976401E-001, -3.1417838584287823E-002,
      2.0528011546134595E-002, -3.3630346024393937E-002, 4.6679798741108718E-002,
      1.7530240854979925E-001, 1.7844697870394516E-001, 1.3828647182431814E-001,
      1.6240726861899379E-001, 1.7898888494748660E-001, 1.4616813940625367E-001,
      1.9950670996141148E-001, 3.3343205576436380E-001, 1.5821214008611131E-001,
      -3.8095185537686921E-001, -2.0184839219705636E-001,
      1.3129407821216750E+000, 1.5097061090315063E+000, -2.5989851697751520E+000,
      -5.3608209887849982E+000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 1.1141460807998201E-003, 2.9893265713566539E-003,
      5.4895923536016841E-003, 8.7378355279600032E-003, 1.2997511828522074E-002,
      1.8633603271577138E-002, 2.6132701485727394E-002, 3.6135981626736598E-002,
      4.9382936813032463E-002, 6.6686176140279288E-002, 8.8733660255192193E-002,
      1.1571165916697089E-001, 1.4640984443289065E-001, 1.7645111710895783E-001,
      1.9471802687302331E-001, 1.7647989007946679E-001, 1.4647104280725104E-001,
      1.1581216380301712E-001, 8.8890066548168750E-002, 6.6910688839154450E-002,
      4.9705881022450292E-002, 3.6580114183300018E-002, 2.6748172235403912E-002,
      1.9462670387787710E-002, 1.4109691154816988E-002, 1.0199423976007888E-002,
      7.3707390106315480E-003, 5.3115014570213094E-003, 3.8130405285331050E-003,
      2.7526316619487573E-003, 1.9778996049961890E-003, 1.4176920470206778E-003,
      1.0193092044854062E-003, 7.3404062497553957E-004, 5.2711429857755707E-004,
      3.7735228904067991E-004, 2.7293872340587879E-004, 1.8988074429557707E-004,
      1.3593101394017193E-004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 4.4565843231992803E-003, 7.5007219622273370E-003,
      1.0001063128980119E-002, 1.2992972697433275E-002, 1.7038705202248282E-002,
      2.2544365772220253E-002, 2.9996392856601019E-002, 4.0013120564036816E-002,
      5.2987820745183462E-002, 6.9212957308987327E-002, 8.8189936459651647E-002,
      1.0791199564711482E-001, 1.2279274106367902E-001, 1.2016509070426867E-001,
      7.3067639056261915E-002, -7.2952547174226112E-002,
      -1.2003538908886301E-001, -1.2263551601693570E-001,
      -1.0768838901939354E-001, -8.7917510836057172E-002,
      -6.8819231266816619E-002, -5.2503067356601110E-002,
      -3.9327767791584423E-002, -2.9142007390464811E-002,
      -2.1411916931882886E-002, -1.5641068715236402E-002,
      -1.1314739861505357E-002, -8.2369502144409544E-003,
      -5.9938437139528174E-003, -4.2416354663373903E-003,
      -3.0989282278102732E-003, -2.2408302319020445E-003,
      -1.5935313701410869E-003, -1.1410743180394665E-003,
      -8.2770530559193021E-004, -5.9904803814750873E-004,
      -4.1765426253920457E-004, -3.3223191644120683E-004,
      -2.1579892142162055E-004, -1.5424054381268617E-004, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 } ;

    helikopterExercise4_last_DWork.x_star_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise4_last_DWork.x_star_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise4_last_DWork.x_star_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopterExercise4_last_terminate();
}

RT_MODEL_helikopterExercise4_la *helikopterExercise4_last(void)
{
  helikopterExercise4_last_initialize(1);
  return helikopterExercise4_last_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
