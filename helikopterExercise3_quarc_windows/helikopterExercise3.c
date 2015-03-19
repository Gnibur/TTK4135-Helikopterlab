/*
 * helikopterExercise3.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterExercise3.mdl".
 *
 * Model version              : 1.64
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Fri Mar 13 12:17:06 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopterExercise3.h"
#include "helikopterExercise3_private.h"
#include <stdio.h>
#include "helikopterExercise3_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopterExercise3 helikopterExercise3_B;

/* Continuous states */
ContinuousStates_helikopterExer helikopterExercise3_X;

/* Block states (auto storage) */
D_Work_helikopterExercise3 helikopterExercise3_DWork;

/* Real-time model */
RT_MODEL_helikopterExercise3 helikopterExercise3_M_;
RT_MODEL_helikopterExercise3 *helikopterExercise3_M = &helikopterExercise3_M_;

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
  helikopterExercise3_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopterExercise3_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Sum1_p[4];
  real_T rtb_Saturation_c;
  real_T rtb_Gain2;
  real_T rtb_Gain1;
  real_T rtb_degrad[6];
  real_T tmp[6];
  int32_T tmp_0;
  int32_T tmp_1;
  real_T tmp_2;
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
    /* set solver stop time */
    if (!(helikopterExercise3_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopterExercise3_M->solverInfo,
                            ((helikopterExercise3_M->Timing.clockTickH0 + 1) *
        helikopterExercise3_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopterExercise3_M->solverInfo,
                            ((helikopterExercise3_M->Timing.clockTick0 + 1) *
        helikopterExercise3_M->Timing.stepSize0 +
        helikopterExercise3_M->Timing.clockTickH0 *
        helikopterExercise3_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopterExercise3_M)) {
    helikopterExercise3_M->Timing.t[0] = rtsiGetT
      (&helikopterExercise3_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise3_B.VandringLavpass =
    helikopterExercise3_P.VandringLavpass_C*
    helikopterExercise3_X.VandringLavpass_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
      if (!(++helikopterExercise3_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopterExercise3_DWork.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise3_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise3_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopterExercise3_M->Timing.t[1];
          u[1] = helikopterExercise3_B.VandringLavpass;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise3_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopterExercise3_DWork.ToFile_IWORK.Count)*2)+1 >=
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

    /* S-Function Block: helikopterExercise3/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder
        (helikopterExercise3_DWork.HILInitialize_Card,
         helikopterExercise3_P.HILReadEncoder_Channels, 3,
         &helikopterExercise3_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 = helikopterExercise3_DWork.HILReadEncoder_Buffer
          [0];
        rtb_HILReadEncoder_o2 = helikopterExercise3_DWork.HILReadEncoder_Buffer
          [1];
        rtb_HILReadEncoder_o3 = helikopterExercise3_DWork.HILReadEncoder_Buffer
          [2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopterExercise3_B.KalibrerElev = helikopterExercise3_P.KalibrerElev_Gain
      * rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopterExercise3_B.Add = helikopterExercise3_B.KalibrerElev +
      helikopterExercise3_P.Constant_Value;

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterExercise3_B.KalibrerPitch =
      helikopterExercise3_P.KalibrerPitch_Gain * rtb_HILReadEncoder_o2;
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopterExercise3_X.Integrator_CSTATE >=
      helikopterExercise3_P.Integrator_UpperSat ) {
    helikopterExercise3_X.Integrator_CSTATE =
      helikopterExercise3_P.Integrator_UpperSat;
  } else if (helikopterExercise3_X.Integrator_CSTATE <=
             helikopterExercise3_P.Integrator_LowerSat ) {
    helikopterExercise3_X.Integrator_CSTATE =
      helikopterExercise3_P.Integrator_LowerSat;
  }

  rtb_Saturation_c = helikopterExercise3_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopterExercise3_B.KalibrerVandring =
      helikopterExercise3_P.KalibrerVandring_Gain * rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopterExercise3_P.VandringDeriv_D*
    helikopterExercise3_B.KalibrerVandring;
  rtb_VandringDeriv += helikopterExercise3_P.VandringDeriv_C*
    helikopterExercise3_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopterExercise3_P.TransferFcn4_D*
    helikopterExercise3_B.KalibrerPitch;
  rtb_Gain2 += helikopterExercise3_P.TransferFcn4_C*
    helikopterExercise3_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopterExercise3_P.TransferFcn5_D*
    helikopterExercise3_B.KalibrerElev;
  rtb_Gain1 += helikopterExercise3_P.TransferFcn5_C*
    helikopterExercise3_X.TransferFcn5_CSTATE;

  /* Gain: '<Root>/deg->rad' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtdeg->radInport1'
   */
  tmp[0] = helikopterExercise3_B.VandringLavpass;
  tmp[1] = rtb_VandringDeriv;
  tmp[2] = helikopterExercise3_B.KalibrerPitch;
  tmp[3] = rtb_Gain2;
  tmp[4] = helikopterExercise3_B.Add;
  tmp[5] = rtb_Gain1;
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_degrad[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_degrad[tmp_0] += helikopterExercise3_P.degrad_Gain[6 * tmp_1 + tmp_0] *
        tmp[tmp_1];
    }
  }

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<Root>/elevation'
   */
  rtb_Gain1 = helikopterExercise3_P.elevation_Value - rtb_degrad[4];

  /* Gain: '<S1>/K_ei' */
  helikopterExercise3_B.K_ei = helikopterExercise3_P.K_ei_Gain * rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ed'
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_c = (helikopterExercise3_P.K_ep_Gain * rtb_Gain1 +
                      rtb_Saturation_c) + helikopterExercise3_P.K_ed_Gain *
    rtb_degrad[5];

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_c = rt_SATURATE(rtb_Saturation_c,
    helikopterExercise3_P.Saturation_LowerSat,
    helikopterExercise3_P.Saturation_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
  }

  /* FromWorkspace: '<Root>/u_star' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise3_DWork.u_star_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise3_DWork.u_star_PWORK.TimePtr;
    int_T currTimeIndex = helikopterExercise3_DWork.u_star_IWORK.PrevIndex;
    real_T t = helikopterExercise3_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helikopterExercise3_DWork.u_star_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Gain1 = pDataValues[currTimeIndex];
        } else {
          rtb_Gain1 = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Gain1 = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* FromWorkspace: '<Root>/x_star' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise3_DWork.x_star_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise3_DWork.x_star_PWORK.TimePtr;
    int_T currTimeIndex = helikopterExercise3_DWork.x_star_IWORK.PrevIndex;
    real_T t = helikopterExercise3_M->Timing.t[0];

    /* get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helikopterExercise3_DWork.x_star_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Sum1_p[0] = pDataValues[currTimeIndex];
          pDataValues += 141;
          rtb_Sum1_p[1] = pDataValues[currTimeIndex];
          pDataValues += 141;
          rtb_Sum1_p[2] = pDataValues[currTimeIndex];
          pDataValues += 141;
          rtb_Sum1_p[3] = pDataValues[currTimeIndex];
          pDataValues += 141;
        } else {
          rtb_Sum1_p[0] = pDataValues[currTimeIndex + 1];
          pDataValues += 141;
          rtb_Sum1_p[1] = pDataValues[currTimeIndex + 1];
          pDataValues += 141;
          rtb_Sum1_p[2] = pDataValues[currTimeIndex + 1];
          pDataValues += 141;
          rtb_Sum1_p[3] = pDataValues[currTimeIndex + 1];
          pDataValues += 141;
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum1_p[0] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum1_p[1] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum1_p[2] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum1_p[3] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* Gain: '<S3>/Gain' incorporates:
   *  Sum: '<S3>/Sum1'
   */
  rtb_Sum1_p[0] = rtb_degrad[0] - rtb_Sum1_p[0];
  tmp_2 = helikopterExercise3_P.Gain_Gain[0] * rtb_Sum1_p[0];
  rtb_Sum1_p[1] = rtb_degrad[1] - rtb_Sum1_p[1];
  tmp_2 += helikopterExercise3_P.Gain_Gain[1] * rtb_Sum1_p[1];
  rtb_Sum1_p[2] = rtb_degrad[3] - rtb_Sum1_p[2];
  tmp_2 += helikopterExercise3_P.Gain_Gain[2] * rtb_Sum1_p[2];
  rtb_Sum1_p[3] = rtb_degrad[4] - rtb_Sum1_p[3];
  rtb_Gain2 = helikopterExercise3_P.Gain_Gain[3] * rtb_Sum1_p[3] + tmp_2;

  /* Sum: '<S3>/Sum' */
  rtb_Gain1 -= rtb_Gain2;

  /* Sum: '<S4>/Sum' incorporates:
   *  Gain: '<S4>/K_pd'
   *  Gain: '<S4>/K_pp'
   *  Saturate: '<S4>/Saturation'
   *  Sum: '<S4>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Gain1,
    helikopterExercise3_P.Saturation_LowerSat_f,
    helikopterExercise3_P.Saturation_UpperSat_d) - rtb_degrad[2]) *
    helikopterExercise3_P.K_pp_Gain - helikopterExercise3_P.K_pd_Gain *
    rtb_degrad[3];

  /* Gain: '<S5>/Gain2' incorporates:
   *  Sum: '<S5>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_c - rtb_Gain1) * helikopterExercise3_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopterExercise3_B.SatB = rt_SATURATE(rtb_Gain2,
    helikopterExercise3_P.SatB_LowerSat, helikopterExercise3_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
  }

  /* Gain: '<S5>/Gain1' incorporates:
   *  Sum: '<S5>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_c) * helikopterExercise3_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopterExercise3_B.Sat = rt_SATURATE(rtb_Gain1,
    helikopterExercise3_P.Sat_LowerSat, helikopterExercise3_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopterExercise3/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopterExercise3_DWork.HILWriteAnalog_Buffer[0] =
        helikopterExercise3_B.SatB;
      helikopterExercise3_DWork.HILWriteAnalog_Buffer[1] =
        helikopterExercise3_B.Sat;
      result = hil_write_analog(helikopterExercise3_DWork.HILInitialize_Card,
        helikopterExercise3_P.HILWriteAnalog_Channels, 2,
        &helikopterExercise3_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopterExercise3_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
    rt_ertODEUpdateContinuousStates(&helikopterExercise3_M->solverInfo);
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
  if (!(++helikopterExercise3_M->Timing.clockTick0)) {
    ++helikopterExercise3_M->Timing.clockTickH0;
  }

  helikopterExercise3_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopterExercise3_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopterExercise3_M->Timing.clockTick1)) {
      ++helikopterExercise3_M->Timing.clockTickH1;
    }

    helikopterExercise3_M->Timing.t[1] =
      helikopterExercise3_M->Timing.clockTick1 *
      helikopterExercise3_M->Timing.stepSize1 +
      helikopterExercise3_M->Timing.clockTickH1 *
      helikopterExercise3_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopterExercise3_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->VandringLavpass_CSTATE = helikopterExercise3_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->VandringLavpass_CSTATE += (helikopterExercise3_P.VandringLavpass_A)*
      helikopterExercise3_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopterExercise3_X.Integrator_CSTATE <=
            helikopterExercise3_P.Integrator_LowerSat );
    usat = ( helikopterExercise3_X.Integrator_CSTATE >=
            helikopterExercise3_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopterExercise3_B.K_ei > 0)) ||
        (usat && (helikopterExercise3_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopterExer *)
        helikopterExercise3_M->ModelData.derivs)->Integrator_CSTATE =
        helikopterExercise3_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopterExer *)
        helikopterExercise3_M->ModelData.derivs)->Integrator_CSTATE = 0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->VandringDeriv_CSTATE = helikopterExercise3_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->VandringDeriv_CSTATE += (helikopterExercise3_P.VandringDeriv_A)*
      helikopterExercise3_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->TransferFcn4_CSTATE = helikopterExercise3_B.KalibrerPitch;
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->TransferFcn4_CSTATE += (helikopterExercise3_P.TransferFcn4_A)*
      helikopterExercise3_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->TransferFcn5_CSTATE = helikopterExercise3_B.KalibrerElev;
    ((StateDerivatives_helikopterExer *) helikopterExercise3_M->ModelData.derivs)
      ->TransferFcn5_CSTATE += (helikopterExercise3_P.TransferFcn5_A)*
      helikopterExercise3_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopterExercise3_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopterExercise3_P.Integrator_UpperSat = rtInf;
  helikopterExercise3_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopterExercise3_M, 0,
                sizeof(RT_MODEL_helikopterExercise3));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopterExercise3_M->solverInfo,
                          &helikopterExercise3_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopterExercise3_M->solverInfo, &rtmGetTPtr
                (helikopterExercise3_M));
    rtsiSetStepSizePtr(&helikopterExercise3_M->solverInfo,
                       &helikopterExercise3_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopterExercise3_M->solverInfo,
                 &helikopterExercise3_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopterExercise3_M->solverInfo,
                         &helikopterExercise3_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopterExercise3_M->solverInfo,
      &helikopterExercise3_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopterExercise3_M->solverInfo,
                          (&rtmGetErrorStatus(helikopterExercise3_M)));
    rtsiSetRTModelPtr(&helikopterExercise3_M->solverInfo, helikopterExercise3_M);
  }

  rtsiSetSimTimeStep(&helikopterExercise3_M->solverInfo, MAJOR_TIME_STEP);
  helikopterExercise3_M->ModelData.intgData.f[0] =
    helikopterExercise3_M->ModelData.odeF[0];
  helikopterExercise3_M->ModelData.contStates = ((real_T *)
    &helikopterExercise3_X);
  rtsiSetSolverData(&helikopterExercise3_M->solverInfo, (void *)
                    &helikopterExercise3_M->ModelData.intgData);
  rtsiSetSolverName(&helikopterExercise3_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopterExercise3_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopterExercise3_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopterExercise3_M->Timing.sampleTimes =
      (&helikopterExercise3_M->Timing.sampleTimesArray[0]);
    helikopterExercise3_M->Timing.offsetTimes =
      (&helikopterExercise3_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopterExercise3_M->Timing.sampleTimes[0] = (0.0);
    helikopterExercise3_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopterExercise3_M->Timing.offsetTimes[0] = (0.0);
    helikopterExercise3_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopterExercise3_M, &helikopterExercise3_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopterExercise3_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopterExercise3_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopterExercise3_M, -1);
  helikopterExercise3_M->Timing.stepSize0 = 0.001;
  helikopterExercise3_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopterExercise3_M->Sizes.checksums[0] = (2016983662U);
  helikopterExercise3_M->Sizes.checksums[1] = (3794230979U);
  helikopterExercise3_M->Sizes.checksums[2] = (632481629U);
  helikopterExercise3_M->Sizes.checksums[3] = (106493464U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopterExercise3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopterExercise3_M->extModeInfo,
      &helikopterExercise3_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopterExercise3_M->extModeInfo,
                        helikopterExercise3_M->Sizes.checksums);
    rteiSetTPtr(helikopterExercise3_M->extModeInfo, rtmGetTPtr
                (helikopterExercise3_M));
  }

  helikopterExercise3_M->solverInfoPtr = (&helikopterExercise3_M->solverInfo);
  helikopterExercise3_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopterExercise3_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopterExercise3_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopterExercise3_M->ModelData.blockIO = ((void *) &helikopterExercise3_B);

  {
    helikopterExercise3_B.VandringLavpass = 0.0;
    helikopterExercise3_B.KalibrerElev = 0.0;
    helikopterExercise3_B.Add = 0.0;
    helikopterExercise3_B.KalibrerPitch = 0.0;
    helikopterExercise3_B.KalibrerVandring = 0.0;
    helikopterExercise3_B.K_ei = 0.0;
    helikopterExercise3_B.SatB = 0.0;
    helikopterExercise3_B.Sat = 0.0;
  }

  /* parameters */
  helikopterExercise3_M->ModelData.defaultParam = ((real_T *)
    &helikopterExercise3_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopterExercise3_X;
    helikopterExercise3_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopterExercise3_X, 0,
                  sizeof(ContinuousStates_helikopterExer));
  }

  /* states (dwork) */
  helikopterExercise3_M->Work.dwork = ((void *) &helikopterExercise3_DWork);
  (void) memset((void *)&helikopterExercise3_DWork, 0,
                sizeof(D_Work_helikopterExercise3));

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise3_DWork.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise3_DWork.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  helikopterExercise3_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopterExercise3_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopterExercise3_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopterExercise3_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopterExercise3_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopterExercise3_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise3/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopterExercise3_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopterExercise3_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopterExercise3_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopterExercise3_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopterExercise3_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopterExercise3_P.HILInitialize_AOExit && is_switching)) {
      helikopterExercise3_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise3_P.HILInitialize_AOFinal;
      helikopterExercise3_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise3_P.HILInitialize_AOFinal;
      helikopterExercise3_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise3_P.HILInitialize_AOFinal;
      helikopterExercise3_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise3_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopterExercise3_DWork.HILInitialize_Card,
        &helikopterExercise3_P.HILInitialize_AOChannels[0], 4U,
        &helikopterExercise3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
      }
    }

    hil_close(helikopterExercise3_DWork.HILInitialize_Card);
    helikopterExercise3_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopterExercise3_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise3_DWork.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error writing header for ans to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      helikopterExercise3_DWork.ToFile_PWORK.FilePtr = (NULL);
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
  helikopterExercise3_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopterExercise3_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopterExercise3_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopterExercise3_M->Sizes.numY = (0);/* Number of model outputs */
  helikopterExercise3_M->Sizes.numU = (0);/* Number of model inputs */
  helikopterExercise3_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopterExercise3_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopterExercise3_M->Sizes.numBlocks = (48);/* Number of blocks */
  helikopterExercise3_M->Sizes.numBlockIO = (8);/* Number of block outputs */
  helikopterExercise3_M->Sizes.numBlockPrms = (153);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise3_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopterExercise3_X.Integrator_CSTATE = helikopterExercise3_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopterExercise3_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopterExercise3_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopterExercise3_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise3/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("sensoray_model_626", "0",
                      &helikopterExercise3_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopterExercise3_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopterExercise3_P.HILInitialize_CKPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_ClockModes =
          &helikopterExercise3_DWork.HILInitialize_ClockModes[0];
        for (i1=0; i1 < 6; i1++) {
          dw_ClockModes[i1] = helikopterExercise3_P.HILInitialize_CKModes;
        }
      }

      result = hil_set_clock_mode(helikopterExercise3_DWork.HILInitialize_Card,
                                  (t_clock *)
        &helikopterExercise3_P.HILInitialize_CKChannels[0], 6U, (t_clock_mode *)
        &helikopterExercise3_DWork.HILInitialize_ClockModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopterExercise3_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
      return;
    }

    if ((helikopterExercise3_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopterExercise3_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helikopterExercise3_DWork.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMinimums[i1] = helikopterExercise3_P.HILInitialize_AILow;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helikopterExercise3_DWork.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMaximums[i1] = helikopterExercise3_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helikopterExercise3_DWork.HILInitialize_Card,
         &helikopterExercise3_P.HILInitialize_AIChannels[0], 16U,
         &helikopterExercise3_DWork.HILInitialize_AIMinimums[0],
         &helikopterExercise3_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise3_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopterExercise3_P.HILInitialize_AOPEnter && is_switching)) {
      helikopterExercise3_DWork.HILInitialize_AOMinimums[0] =
        helikopterExercise3_P.HILInitialize_AOLow;
      helikopterExercise3_DWork.HILInitialize_AOMinimums[1] =
        helikopterExercise3_P.HILInitialize_AOLow;
      helikopterExercise3_DWork.HILInitialize_AOMinimums[2] =
        helikopterExercise3_P.HILInitialize_AOLow;
      helikopterExercise3_DWork.HILInitialize_AOMinimums[3] =
        helikopterExercise3_P.HILInitialize_AOLow;
      helikopterExercise3_DWork.HILInitialize_AOMaximums[0] =
        helikopterExercise3_P.HILInitialize_AOHigh;
      helikopterExercise3_DWork.HILInitialize_AOMaximums[1] =
        helikopterExercise3_P.HILInitialize_AOHigh;
      helikopterExercise3_DWork.HILInitialize_AOMaximums[2] =
        helikopterExercise3_P.HILInitialize_AOHigh;
      helikopterExercise3_DWork.HILInitialize_AOMaximums[3] =
        helikopterExercise3_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopterExercise3_DWork.HILInitialize_Card,
         &helikopterExercise3_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise3_DWork.HILInitialize_AOMinimums[0],
         &helikopterExercise3_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise3_P.HILInitialize_AOStart && !is_switching) ||
        (helikopterExercise3_P.HILInitialize_AOEnter && is_switching)) {
      helikopterExercise3_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise3_P.HILInitialize_AOInitial;
      helikopterExercise3_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise3_P.HILInitialize_AOInitial;
      helikopterExercise3_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise3_P.HILInitialize_AOInitial;
      helikopterExercise3_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise3_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopterExercise3_DWork.HILInitialize_Card,
        &helikopterExercise3_P.HILInitialize_AOChannels[0], 4U,
        &helikopterExercise3_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise3_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise3_M,
                        "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise3_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopterExercise3_DWork.ToFile_IWORK.Count = 0;
    helikopterExercise3_DWork.ToFile_IWORK.Decimation = -1;
    helikopterExercise3_DWork.ToFile_PWORK.FilePtr = fp;
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      2.9253550210096735E-001, 3.6566937762617234E-002, 5.2359877559829837E-001,
      5.2359877559829726E-001, 5.2359877559829793E-001, 5.2359877559829826E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829393E-001,
      5.2359877559829482E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829759E-001, 5.2359877559829526E-001,
      5.2359877559829804E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 4.8525821529535401E-001,
      3.9579764804160383E-001, 3.0767547240042847E-001, 1.4396515428661602E-001,
      -3.8868990282269977E-002, -1.6638694992517056E-001,
      -2.3482225198551832E-001, -2.8833194095600462E-001,
      -3.5378622094461132E-001, -4.2015416237802045E-001,
      -4.6682735553366617E-001, -4.9087711533085931E-001,
      -5.0436066749931974E-001, -5.1653003628552929E-001,
      -5.2359877559829882E-001, -5.2216894584605200E-001,
      -5.1711677207739570E-001, -5.0988203136472443E-001,
      -4.9692489887036201E-001, -4.7846170170730429E-001,
      -4.5742120266403391E-001, -4.3583815656684366E-001,
      -4.1358745229039307E-001, -3.8988763785556618E-001,
      -3.6485417153474714E-001, -3.3941568272474854E-001,
      -3.1434027265584580E-001, -2.8977739198258901E-001,
      -2.6560516958605385E-001, -2.4189131862219479E-001,
      -2.1892992669283945E-001, -1.9698709503139403E-001,
      -1.7614821149016649E-001, -1.5639338093144722E-001,
      -1.3773157936221789E-001, -1.2022849827480481E-001,
      -1.0394311100396064E-001, -8.8879971844661787E-002,
      -7.5004447135339367E-002, -6.2281409861114738E-002,
      -5.0688313887411157E-002, -4.0200563151115004E-002,
      -3.0777559171217831E-002, -2.2365592481432352E-002,
      -1.4909032770836656E-002, -8.3556352689666234E-003,
      -2.6536559286987633E-003, 2.2519249844980180E-003, 6.4188640470197724E-003,
      9.9044394726117771E-003, 1.2763502601698879E-002, 1.5048846310942859E-002,
      1.6812103799925072E-002, 1.8103683683243364E-002, 1.8972033036834628E-002,
      1.9463144188227234E-002, 1.9620475076426152E-002, 1.9484916045897343E-002,
      1.9094805811813997E-002, 1.8486312446061613E-002, 1.7693920221605031E-002,
      1.6750160044295175E-002, 1.5684512888115662E-002, 1.4523050211504546E-002,
      1.3290237701381834E-002, 1.2011393490285386E-002, 1.0712054473348932E-002,
      9.4135435749134605E-003, 8.1304126090638983E-003, 6.8758916682217239E-003,
      5.6712162929760217E-003, 4.5449906165087647E-003, 3.5173612993384812E-003,
      2.5874917078657299E-003, 1.7492101661484114E-003, 1.0257958473993995E-003,
      4.7457164521870162E-004, 1.3318906902382239E-004, -1.2137712376410970E-016,
      3.2607753421329408E-016, 3.2607753421329408E-016, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    } ;

    helikopterExercise3_DWork.u_star_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise3_DWork.u_star_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise3_DWork.u_star_IWORK.PrevIndex = 0;
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues[] = { 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897922E+000,
      3.1415926535897900E+000, 3.1415926535897927E+000, 3.1358533890859204E+000,
      3.1243748600781736E+000, 3.1071570665665509E+000, 3.0842000085510550E+000,
      3.0555036860316851E+000, 3.0210680990084429E+000, 2.9808932474813266E+000,
      2.9349791314503344E+000, 2.8833257509154677E+000, 2.8259331058767287E+000,
      2.7628011963341090E+000, 2.6939300222876210E+000, 2.6193195837372620E+000,
      2.5389698806830276E+000, 2.4528809131249174E+000, 2.3610526810629331E+000,
      2.2634851844970716E+000, 2.1601784234273387E+000, 2.0511323978537277E+000,
      1.9363471077762475E+000, 1.8165747580220808E+000, 1.6942286578913686E+000,
      1.5718970063776350E+000, 1.4510384843298092E+000, 1.3322074676957354E+000,
      1.2159623095607961E+000, 1.1032112429392806E+000, 9.9486919825164122E-001,
      8.9150880499628016E-001, 7.9340995425070748E-001, 7.0077988349758824E-001,
      6.1380577804391301E-001, 5.3257093709410674E-001, 4.5707536064818260E-001,
      3.8731904870612632E-001, 3.2327394938677245E-001, 2.6481639851579719E-001,
      2.1174295387144809E-001, 1.6382728109960165E-001, 1.2084781115248626E-001,
      8.2574692735378916E-002, 4.8752249750204625E-002, 1.9102899401296527E-002,
      -6.6555448251006573E-003, -2.8798990987201314E-002,
      -4.7598135080634242E-002, -6.3322442066317611E-002,
      -7.6238921792756972E-002, -8.6606821051787075E-002,
      -9.4673688583468960E-002, -1.0067499924564967E-001,
      -1.0483515574458237E-001, -1.0736731562843448E-001,
      -1.0847196248475439E-001, -1.0833576766993153E-001,
      -1.0713168428070043E-001, -1.0501976744491905E-001,
      -1.0214764627532906E-001, -9.8650250635573741E-002,
      -9.4649475479723452E-002, -9.0254784380256428E-002,
      -8.5564778995128196E-002, -8.0668364062470979E-002,
      -7.5644207549631179E-002, -7.0559420138837495E-002,
      -6.5470323547670323E-002, -6.0426336454573203E-002,
      -5.5473174968744175E-002, -5.0650316631950451E-002,
      -4.5984576816235620E-002, -4.1489786294499856E-002,
      -3.7178239284569231E-002, -3.3072260282485985E-002,
      -2.9196254573170606E-002, -2.5551445707358150E-002,
      -2.2107907937982167E-002, -1.8840505871074636E-002,
      -1.5774408799594689E-002, -1.2965833281753433E-002,
      -1.0410127302073425E-002, -7.9940814493383791E-003,
      -5.6078667972605332E-003, -3.3229543432832880E-003,
      -1.3619446621948980E-003, 2.2356282944416942E-004, 1.7841080247865825E-003,
      3.7378111222137167E-003, 5.8924015645241933E-003, 7.4087031239602336E-003,
      7.8971534291596198E-003, 8.4376342992754701E-003, 1.0645706171359865E-002,
      1.4157289116225682E-002, 1.5982367864786993E-002, 1.4127268913963667E-002,
      1.1839689540199634E-002, 1.5291374670310380E-002, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.4884255795919775E-016,
      -8.4545807725769899E-016, -2.2957058015496235E-002,
      -4.5914116030992207E-002, -6.8871174046487987E-002,
      -9.1828232061983955E-002, -1.1478529007748088E-001,
      -1.3774234809297764E-001, -1.6069940610847050E-001,
      -1.8365646412396783E-001, -2.0661352213946446E-001,
      -2.2957058015496137E-001, -2.5252763817045865E-001,
      -2.7548469618595389E-001, -2.9844175420144919E-001,
      -3.2139881221694594E-001, -3.4435587023244280E-001,
      -3.6731292824793993E-001, -3.9026998626343623E-001,
      -4.1322704427893464E-001, -4.3618410229443250E-001,
      -4.5914116030992874E-001, -4.7908939901665770E-001,
      -4.8938440052283777E-001, -4.8932660605494743E-001,
      -4.8343408819130024E-001, -4.7532406653629661E-001,
      -4.6498063253975686E-001, -4.5100426648607533E-001,
      -4.3336817875054912E-001, -4.1344157302145518E-001,
      -3.9239540298228054E-001, -3.7052028301246764E-001,
      -3.4789642181470060E-001, -3.2493936379920502E-001,
      -3.0198230578370933E-001, -2.7902524776821480E-001,
      -2.5618039727740810E-001, -2.3383020348390859E-001,
      -2.1229377857739740E-001, -1.9166269108738235E-001,
      -1.7191787978846765E-001, -1.5309247366842485E-001,
      -1.3528977194069211E-001, -1.1859740139562251E-001,
      -1.0303377690559529E-001, -8.8573784648401488E-002,
      -7.5196576373743301E-002, -6.2897227942743023E-002,
      -5.1665918905756468E-002, -4.1471597036129258E-002,
      -3.2267470126728037E-002, -2.4005242648723193E-002,
      -1.6640625995732992E-002, -1.0128639535411292E-002,
      -4.4185874252750590E-003, 5.4477925929091470E-004, 4.8163335569151892E-003,
      8.4476673431339386E-003, 1.1488484678360061E-002, 1.3989582559020020E-002,
      1.6003100623399994E-002, 1.7578764397864193E-002, 1.8760021540513129E-002,
      1.9585659730627529E-002, 2.0096626051363808E-002, 2.0339149643169050E-002,
      2.0356386364668771E-002, 2.0175948372383790E-002, 1.9812645943312456E-002,
      1.9291433347172572E-002, 1.8662959262859186E-002, 1.7979162086941336E-002,
      1.7246188039723668E-002, 1.6423916008331870E-002, 1.5504022837259573E-002,
      1.4579235463248120E-002, 1.3774151077504798E-002, 1.3069608267628381E-002,
      1.2264388285920681E-002, 1.1234302071364515E-002, 1.0222823918717368E-002,
      9.6641834109388435E-003, 9.5448586083130040E-003, 9.1396498159069424E-003,
      7.8440387243548168E-003, 6.3420299665552368E-003, 6.2421807813712055E-003,
      7.8148123897087887E-003, 8.6183617692433493E-003, 6.0652062377445975E-003,
      1.9538012207994640E-003, 2.1619234804629160E-003, 8.8322874883359909E-003,
      1.4046331779460056E-002, 7.3003149942470454E-003, -7.4203958032933613E-003,
      -9.1503174950535952E-003, 1.3806740520442914E-002, 3.6763798535939189E-002,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      4.8977584246598528E-018, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829715E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829604E-001,
      5.2359877559829626E-001, 5.2359877559829560E-001, 5.2359877559829882E-001,
      5.2359877559829870E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829537E-001, 5.2359877559829882E-001, 5.2359877559829604E-001,
      5.2359877559829349E-001, 5.2359877559829071E-001, 5.2359877559829882E-001,
      4.5497438544326818E-001, 2.3480579174300403E-001, -1.3181616130222072E-003,
      -1.3439505778634472E-001, -1.8497132366739633E-001,
      -2.3591042773923526E-001, -3.1876942368142758E-001,
      -4.0223943061136952E-001, -4.5448098595821079E-001,
      -4.8001572571304862E-001, -4.9892220616985139E-001,
      -5.1599930681277228E-001, -5.2359877559829793E-001,
      -5.2359877559829882E-001, -5.2359877559829882E-001,
      -5.2103957474160922E-001, -5.0975756984026288E-001,
      -4.9119733478926075E-001, -4.7054862786601820E-001,
      -4.5033466455293886E-001, -4.2936510366184844E-001,
      -4.0603952042510771E-001, -3.8071536750618362E-001,
      -3.5497121283349131E-001, -3.2979984786883065E-001,
      -3.0510398452652376E-001, -2.8052042969879565E-001,
      -2.5616085720400816E-001, -2.3250951604462289E-001,
      -2.0992539971629917E-001, -1.8844279560026839E-001,
      -1.6797031482227820E-001, -1.4852374093559817E-001,
      -1.3023342500818347E-001, -1.1320321258636004E-001,
      -9.7424530557390837E-002, -8.2822542983050376E-002,
      -6.9354193053303440E-002, -5.7044408176380831E-002,
      -4.5923810988440325E-002, -3.5937341035027749E-002,
      -2.6941814283865679E-002, -1.8830947116111194E-002,
      -1.1653990669397391E-002, -5.5314167711429038E-003,
      -3.9313078645259489E-004, 4.1153841127192434E-003, 8.2861099582180932E-003,
      1.1887685128518248E-002, 1.4334078034717939E-002, 1.5595873122186677E-002,
      1.6717486770712379E-002, 1.8754172618928359E-002, 2.0980690893831748E-002,
      2.1092316637192614E-002, 1.8362161142074274E-002, 1.6069034296949534E-002,
      1.8365253780539919E-002, 2.3493945972434729E-002, 2.3069538000580897E-002,
      1.2741331475281918E-002, 2.7215299325875717E-003, 9.2418996990945301E-003,
      2.9549970241410463E-002, 3.4257436035187334E-002, 2.2773349734773437E-003,
      -3.5868184156560873E-002, -1.8327151108524807E-002,
      5.8231725916374454E-002, 9.3771886249878922E-002, -4.7467998844156452E-003,
      -1.5213597600178430E-001, -1.1892060406455289E-001,
      1.5386144542208977E-001, 3.3574625042660478E-001, 3.9455616615837220E-002,
      -5.2359877559829882E-001, -5.2359877559829882E-001,
      4.1357184416948251E-001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 2.0943951023931957E+000, -3.2265225704909477E-015,
      3.9211564422617098E-015, -5.6639745197151216E-016,
      -8.3958121246599677E-015, 6.3990163019349049E-015, 3.7677573255933283E-015,
      2.5784521538926860E-016, 6.7380111633214734E-015, -1.1474860746051581E-015,
      1.3169882086829746E-014, 6.1312535756340001E-015, 3.1755489675725616E-016,
      -2.8316912885174545E-015, -2.0605738590171115E-014,
      6.8118839214221906E-015, -1.5061959407786447E-014,
      -9.0043077798903375E-015, -1.3454508952369329E-014,
      3.7034578901097653E-014, -2.7449756062011055E-001,
      -8.8067437480105670E-001, -9.4449581342410438E-001,
      -5.3230758469329142E-001, -2.0230506352421312E-001,
      -2.0375641628735189E-001, -3.3143598376876715E-001,
      -3.3388002771975278E-001, -2.0896622138736781E-001,
      -1.0213895901935201E-001, -7.5625921827210341E-002,
      -6.8308402571700594E-002, -3.0397875142099597E-002,
      -1.9841755022970859E-015, -7.5560877394431530E-015,
      1.0236803426766894E-002, 4.5128019605378829E-002, 7.4240940204008732E-002,
      8.2594827692969358E-002, 8.0855853252304072E-002, 8.3878243564359523E-002,
      9.3302332946962288E-002, 1.0129661167569103E-001, 1.0297661869077210E-001,
      1.0068545985863954E-001, 9.8783453369219246E-002, 9.8334219310909660E-002,
      9.7438289979163928E-002, 9.4605364637552897E-002, 9.0336465313290384E-002,
      8.5930416464118947E-002, 8.1889923111964497E-002, 7.7786295546723255E-002,
      7.3161263709656130E-002, 6.8120849687291596E-002, 6.3114728115867061E-002,
      5.8407950297357986E-002, 5.3873399718992590E-002, 4.9239139507688705E-002,
      4.4482388751760124E-002, 3.9945879813644974E-002, 3.5982107004647171E-002,
      3.2443468671019271E-002, 2.8707825786857984E-002, 2.4490295593020253E-002,
      2.0553143938760454E-002, 1.8034059596691070E-002, 1.6682903381997075E-002,
      1.4406300681199722E-002, 9.7855716247983431E-003, 5.0471803498768723E-003,
      4.4864545941060525E-003, 8.1467433928660987E-003, 8.9060730996133892E-003,
      4.4650297344342980E-004, -1.0920621980474823E-002,
      -9.1725073804990766E-003, 9.1848779343623732E-003, 2.0514768767580918E-002,
      -1.6976318874134264E-003, -4.1312826101191483E-002,
      -4.0079206170776586E-002, 2.6081479066028701E-002, 8.1232282169265188E-002,
      1.8829863175108194E-002, -1.2792040424683770E-001,
      -1.5258207652015177E-001, 7.0164132192144887E-002, 3.0623550809959688E-001,
      1.4216064133401746E-001, -3.9407474453718089E-001,
      -5.8955670446947595E-001, 1.3286148774892692E-001, 1.0911281979465717E+000,
      7.2753922001806082E-001, -1.1851625352430681E+000,
      -2.2522175688565471E+000, -6.7853921765690593E-016,
      3.7486824790711255E+000, 3.2800971691872327E+000, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    } ;

    helikopterExercise3_DWork.x_star_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise3_DWork.x_star_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise3_DWork.x_star_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopterExercise3_terminate();
}

RT_MODEL_helikopterExercise3 *helikopterExercise3(void)
{
  helikopterExercise3_initialize(1);
  return helikopterExercise3_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
