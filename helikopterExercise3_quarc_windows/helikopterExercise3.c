/*
 * helikopterExercise3.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterExercise3.mdl".
 *
 * Model version              : 1.80
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Tue Apr 28 12:33:40 2015
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

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterExercise3_B.KalibrerPitch =
      helikopterExercise3_P.KalibrerPitch_Gain * rtb_HILReadEncoder_o2;

    /* ToFile: '<Root>/To File1' */
    if (rtmIsMajorTimeStep(helikopterExercise3_M)) {
      if (!(++helikopterExercise3_DWork.ToFile1_IWORK.Decimation % 1) &&
          (helikopterExercise3_DWork.ToFile1_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise3_DWork.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise3_DWork.ToFile1_IWORK.Decimation = 0;
          u[0] = helikopterExercise3_M->Timing.t[1];
          u[1] = helikopterExercise3_B.KalibrerPitch;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise3_M,
                              "Error writing to MAT-file pitch.mat");
            return;
          }

          if (((++helikopterExercise3_DWork.ToFile1_IWORK.Count)*2)+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file pitch.mat.\n");
          }
        }
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

  /* Gain: '<S3>/Gain' */

  /* Sum: '<S3>/Sum1' */
  rtb_Sum1_p[0] = rtb_degrad[0] - rtb_Sum1_p[0];
  tmp_2 = helikopterExercise3_P.Gain_Gain[0] * rtb_Sum1_p[0];

  /* Sum: '<S3>/Sum1' */
  rtb_Sum1_p[1] = rtb_degrad[1] - rtb_Sum1_p[1];
  tmp_2 += helikopterExercise3_P.Gain_Gain[1] * rtb_Sum1_p[1];

  /* Sum: '<S3>/Sum1' */
  rtb_Sum1_p[2] = rtb_degrad[2] - rtb_Sum1_p[2];
  tmp_2 += helikopterExercise3_P.Gain_Gain[2] * rtb_Sum1_p[2];

  /* Sum: '<S3>/Sum1' */
  rtb_Sum1_p[3] = rtb_degrad[3] - rtb_Sum1_p[3];
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
  helikopterExercise3_M->Sizes.checksums[0] = (3940357121U);
  helikopterExercise3_M->Sizes.checksums[1] = (2279441428U);
  helikopterExercise3_M->Sizes.checksums[2] = (4036358734U);
  helikopterExercise3_M->Sizes.checksums[3] = (2924317805U);

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
    helikopterExercise3_B.KalibrerPitch = 0.0;
    helikopterExercise3_B.KalibrerElev = 0.0;
    helikopterExercise3_B.Add = 0.0;
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

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) helikopterExercise3_DWork.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "pitch.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error reopening MAT-file pitch.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise3_DWork.ToFile1_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error writing header for ans to MAT-file pitch.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise3_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      helikopterExercise3_DWork.ToFile1_PWORK.FilePtr = (NULL);
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
  helikopterExercise3_M->Sizes.numBlocks = (49);/* Number of blocks */
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

  /* Start for ToFile: '<Root>/To File1' */
  {
    const char *fileName = "pitch.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise3_M,
                        "Error creating .mat file pitch.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise3_M,
                        "Error writing mat file header to file pitch.mat");
      return;
    }

    helikopterExercise3_DWork.ToFile1_IWORK.Count = 0;
    helikopterExercise3_DWork.ToFile1_IWORK.Decimation = -1;
    helikopterExercise3_DWork.ToFile1_PWORK.FilePtr = fp;
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
      1.6455121993178229E-001, 8.2275609965902552E-002, 5.2359877559829371E-001,
      5.2359877559829326E-001, 5.2359877559829859E-001, 5.2359877559829882E-001,
      5.2359877559829870E-001, 5.2359877559829882E-001, 5.2359877559829859E-001,
      5.2359877559829726E-001, 5.2359877559829382E-001, 5.2359877559829815E-001,
      5.2359877559829837E-001, 5.2359877559829882E-001, 4.7672641088458512E-001,
      4.0573593707154687E-001, 2.5668373741646877E-001, 4.1989344833112119E-003,
      -1.9215902420577685E-001, -2.9795864027284069E-001,
      -3.7095644928699822E-001, -4.3928069310872719E-001,
      -4.9095547989911648E-001, -5.1731839163327975E-001,
      -5.2282559327494593E-001, -5.1819765055675293E-001,
      -5.1018177555175037E-001, -4.9124220194208745E-001,
      -4.6319132600830115E-001, -4.3134853989393296E-001,
      -3.9790790860711078E-001, -3.6289861635459164E-001,
      -3.2685695654417118E-001, -2.9098766717717184E-001,
      -2.5622002441247943E-001, -2.2299409035526005E-001,
      -1.9159460673146161E-001, -1.6231389133064281E-001,
      -1.3537762754355351E-001, -1.1088896004955637E-001,
      -8.8862990935106856E-002, -6.9269802331373947E-002,
      -5.2043840128313298E-002, -3.7082331774687484E-002,
      -2.4252716243106848E-002, -1.3404412584526694E-002,
      -4.3766838264000284E-003, 2.9971082619017946E-003, 8.8856317396846449E-003,
      1.3455503164786428E-002, 1.6867984079489893E-002, 1.9276756620555516E-002,
      2.0826364837195521E-002, 2.1651035205005126E-002, 2.1873885260112979E-002,
      2.1606543177297325E-002, 2.0949100290611249E-002, 1.9990301634163547E-002,
      1.8807921170164998E-002, 1.7469289520661926E-002, 1.6031937353419894E-002,
      1.4544315422005141E-002, 1.3046559828020542E-002, 1.1571279084749172E-002,
      1.0144343565654624E-002, 8.7856605318402149E-003, 7.5099212879548745E-003,
      6.3273105434124444E-003, 5.2441708601465142E-003, 4.2636172308684458E-003,
      3.3860987392647546E-003, 2.6099059564871348E-003, 1.9316240750686805E-003,
      1.3465327972450828E-003, 8.4895491925216220E-004, 4.3255642875962791E-004,
      9.0601303033112108E-005, -1.8383613665669108E-004,
      -3.9769830889009253E-004, -5.5780607743870822E-004,
      -6.7074930471780220E-004, -7.4280471850093953E-004,
      -7.7988435327993916E-004, -7.8750916748105033E-004,
      -7.7078571299276266E-004, -7.3437917973171460E-004,
      -6.8252936408872095E-004, -6.1914196259337439E-004,
      -5.4783846050499214E-004, -4.7181200428522078E-004,
      -3.9371938637493503E-004, -3.1612737955520099E-004,
      -2.4214959813691680E-004, -1.7469704914750990E-004,
      -1.1452573617906442E-004, -6.1398159200084929E-005,
      -2.0020037504940381E-005, 2.6584814020533603E-015,
      -1.9565406777024674E-015, -7.9383773642678285E-016,
      -7.9383773642678285E-016, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1415926535897931E+000,
      3.1415926535897913E+000, 3.1415926535897905E+000, 3.1296829880716448E+000,
      3.1058636570353468E+000, 3.0701346604809001E+000, 3.0224959984083042E+000,
      2.9629476708175573E+000, 2.8914896777086629E+000, 2.8081220190816176E+000,
      2.7128446949364253E+000, 2.6056577052730820E+000, 2.4865610500915887E+000,
      2.3555547293919470E+000, 2.2126387431741521E+000, 2.0578130914382151E+000,
      1.8910777741841314E+000, 1.7158252559953666E+000, 1.5388898216098044E+000,
      1.3653779878870778E+000, 1.1977875981466790E+000, 1.0378837259431730E+000,
      8.8750716237956018E-001, 7.4812774676274818E-001, 6.2048851253505621E-001,
      5.0475894382550746E-001, 4.0093904063411206E-001, 3.0897284261696900E-001,
      2.2849741462681819E-001, 1.5890121148830427E-001, 9.9481062778218179E-002,
      4.9495307190841248E-002, 8.1543293650563481E-003, -2.5367877966728913E-002,
      -5.1900440603653106E-002, -7.2251160528682509E-002,
      -8.7196338381118899E-002, -9.7473980700929552E-002,
      -1.0377584926181672E-001, -1.0674109512477566E-001,
      -1.0695324221730093E-001, -1.0493941826019712E-001,
      -1.0117059494684062E-001, -9.6062677633638349E-002,
      -8.9978571384282391E-002, -8.3231036049886187E-002,
      -7.6086011051314198E-002, -6.8766214823655192E-002,
      -6.1454926081218587E-002, -5.4299854413184627E-002,
      -4.7416995018709378E-002, -4.0894382840932789E-002,
      -3.4795690457255168E-002, -2.9163630079930829E-002,
      -2.4023128082447769E-002, -1.9384249777033950E-002,
      -1.5244862647681743E-002, -1.1593034604843377E-002,
      -8.4091696466180380E-003, -5.6678878921713124E-003,
      -3.3396608646860020E-003, -1.3922159384635634E-003,
      2.0827406990313031E-004, 1.4961989288858476E-003, 2.5056987629734511E-003,
      3.2700134729883824E-003, 3.8209757339021382E-003, 4.1886299248609302E-003,
      4.4009605434117468E-003, 4.4837151274602118E-003, 4.4603069689174634E-003,
      4.3517823369976583E-003, 4.1768395572613351E-003, 3.9518939544131906E-003,
      3.6911828768442820E-003, 3.4068911858669993E-003, 3.1092739658982872E-003,
      2.8067922168313005E-003, 2.5063100284794097E-003, 2.2133167742374241E-003,
      1.9320106472970095E-003, 1.6652469786609804E-003, 1.4147820248984659E-003,
      1.1820104242433059E-003, 9.6813437914570252E-004, 7.7269829789025483E-004,
      5.9280195955573468E-004, 4.2700675422587610E-004, 2.7973673014856286E-004,
      1.5330951941931856E-004, 3.2430959319982664E-005, -1.0192258401311583E-004,
      -2.2374682106989996E-004, -2.7467067356281210E-004,
      -3.0086744939804169E-004, -4.8732816090883371E-004,
      -8.0759439143955863E-004, -7.5664963993180739E-004,
      -1.5679180967158093E-004, -2.7660907775340537E-004, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.3475727763221047E-016,
      1.3364922691988183E-016, -4.7638662072596573E-002,
      -9.5277324145193784E-002, -1.4291598621779009E-001,
      -1.9055464829038679E-001, -2.3819331036298308E-001,
      -2.8583197243557873E-001, -3.3347063450817738E-001,
      -3.8110929658077447E-001, -4.2874795865337006E-001,
      -4.7638662072596727E-001, -5.2402528279856364E-001,
      -5.7166394487116101E-001, -6.1930260694375683E-001,
      -6.6694126901635320E-001, -7.0101007275507654E-001,
      -7.0774173754224912E-001, -6.9404733489089099E-001,
      -6.7036155896160565E-001, -6.3961548881401908E-001,
      -6.0150625425444915E-001, -5.5751766246724754E-001,
      -5.1055693691077919E-001, -4.6291827483818077E-001,
      -4.1527961276558356E-001, -3.6786479206856609E-001,
      -3.2190171196061274E-001, -2.7838481255404296E-001,
      -2.3768059484034101E-001, -1.9994302234949649E-001,
      -1.6536391130316802E-001, -1.3408882932714769E-001,
      -1.0613025054769018E-001, -8.1402879700108635E-002,
      -5.9780711409743101E-002, -4.1110569279249519E-002,
      -2.5207474243552342E-002, -1.1860983451847157E-002,
      -8.4858837009217029E-004, 8.0552958284130418E-003, 1.5075293253433177E-002,
      2.0431669252814300E-002, 2.4336424997443453E-002, 2.6990141337596994E-002,
      2.8580099994283244E-002, 2.9279184910628964E-002, 2.9245154969741970E-002,
      2.8620286672147988E-002, 2.7531437577913521E-002, 2.6090448711116848E-002,
      2.4394769534710797E-002, 2.2528241509298424E-002, 2.0562007989928972E-002,
      1.8555513221662056E-002, 1.6557548517408011E-002, 1.4607312171359389E-002,
      1.2735459832899743E-002, 1.0965127017783335E-002, 9.3129081099393978E-003,
      7.7897797048904222E-003, 6.4019600334679315E-003, 5.1516994359341317E-003,
      4.0379993363464988E-003, 3.0572588400580789E-003, 2.2038490436475174E-003,
      1.4706167638339903E-003, 8.4932247420345353E-004, 3.3101833619512397E-004,
      -9.3632634174449812E-005, -4.3409852767641765E-004,
      -6.9977111894733047E-004, -8.9978241139427781E-004,
      -1.0428443102787302E-003, -1.1371667639089387E-003,
      -1.1904688798737204E-003, -1.2099269962667415E-003,
      -1.2019287534065590E-003, -1.1719730169662485E-003,
      -1.1252245077599057E-003, -1.0670546745419180E-003,
      -1.0018598150491680E-003, -9.3108640262143014E-004,
      -8.5550418039057107E-004, -7.8174432502058750E-004,
      -7.1958535333873459E-004, -6.6318082131921997E-004,
      -5.8908009630968382E-004, -5.0570884291799702E-004,
      -4.8351424039810137E-004, -5.3741417333303277E-004,
      -4.8729694822727347E-004, -2.0369540997048580E-004,
      -1.0478710333922169E-004, -7.4584284604324797E-004,
      -1.2810649221225885E-003, 2.0377900602360099E-004, 2.3994313210430329E-003,
      -4.7926907232402393E-004, -7.8073159963509007E-003, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6462595819664960E-017,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829837E-001, 5.2359877559829882E-001, 5.2359877559829859E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829593E-001, 3.7445182437137015E-001,
      7.3988044310112397E-002, -1.5051582367856350E-001,
      -2.6033147733599094E-001, -3.3793150318126897E-001,
      -4.1886038957119731E-001, -4.8348068140478950E-001,
      -5.1614754346173108E-001, -5.2359877559829882E-001,
      -5.2359877559829882E-001, -5.2113852451058229E-001,
      -5.0518237541132860E-001, -4.7829628826248155E-001,
      -4.4738197147733105E-001, -4.1477543429216407E-001,
      -3.8006063599234841E-001, -3.4374589707042719E-001,
      -3.0729405444017255E-001, -2.7177969607603586E-001,
      -2.3765026871165834E-001, -2.0520441033535933E-001,
      -1.7479166556409409E-001, -1.4669190806452931E-001,
      -1.2103775232866404E-001, -9.7863010033771763E-002,
      -7.7157121894867592E-002, -5.8872180554674740E-002,
      -4.2917354055491047E-002, -2.9167121116293916E-002,
      -1.7475310381806337E-002, -7.6836751980776793E-003,
      3.7402468094747168E-004, 6.8679568505112081E-003, 1.1967591610429391E-002,
      1.5837976414110999E-002, 1.8637289586788387E-002, 2.0515097322336751E-002,
      2.1610965096217963E-002, 2.2053478376617832E-002, 2.1959724041828492E-002,
      2.1435139411818637E-002, 2.0573617097492728E-002, 1.9457811241287235E-002,
      1.8159615730784791E-002, 1.6740775942602870E-002, 1.5253591286856020E-002,
      1.3741673035432373E-002, 1.2240730179185118E-002, 1.0779364925371357E-002,
      9.3798672138118988E-003, 8.0589904761502670E-003, 6.8286747608672276E-003,
      5.6967051600948174E-003, 4.6673587894543135E-003, 3.7420766498156635E-003,
      2.9200199469810702E-003, 2.1983335231330267E-003, 1.5724000597712030E-003,
      1.0367025244559685E-003, 5.8584606370275791E-004, 2.1386507251863570E-004,
      -8.7909063481875384E-005, -3.2924490823045088E-004,
      -5.1381506358891336E-004, -6.3934737301849724E-004,
      -7.1655976722493817E-004, -7.7787390492968045E-004,
      -8.3072775966209198E-004, -8.1069803977351488E-004,
      -6.8319218149621443E-004, -6.1994486449837597E-004,
      -8.1444455402151139E-004, -9.1633736751854666E-004,
      -2.4394192024915080E-004, 5.9241669815855952E-004,
      -5.5084077847809777E-004, -3.1170778466182167E-003,
      -1.0871058504873347E-003, 7.0458738211115965E-003, 5.8826510132018390E-003,
      -1.6319989456188101E-002, -2.4132517870295491E-002,
      3.1639931427638532E-002, 8.0542908428039531E-002, -4.4586225381588532E-002,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      2.0943951023931895E+000, 7.2346488159173625E-016, -2.4461319426746430E-015,
      -3.4487305754495795E-015, 1.1176156751535719E-015,
      -3.5465108439355480E-016, 4.3518371215519393E-015,
      -4.3221941493288452E-015, 2.0547038476551659E-015,
      -2.2432101341000978E-015, -3.8929392183791486E-015,
      -2.2159295454315319E-015, -6.6797327362753827E-015,
      -7.3720918322529251E-015, -5.9658780490769736E-001,
      -1.2018551202450227E+000, -8.9801547195470621E-001,
      -4.3926261462972011E-001, -3.1040010338111362E-001,
      -3.2371554555971671E-001, -2.5848116733437509E-001,
      -1.3066744822776952E-001, -2.9804928546273004E-002,
      -1.5522537028070505E-014, 9.8410043508775913E-003, 6.3824596397015515E-002,
      1.0754434859537639E-001, 1.2365726714060765E-001, 1.3042614874066591E-001,
      1.3885919319925943E-001, 1.4525895568768390E-001, 1.4580737052100856E-001,
      1.4205743345653910E-001, 1.3651770945750660E-001, 1.2978343350519797E-001,
      1.2165097908506049E-001, 1.1239902999826933E-001, 1.0261662294346105E-001,
      9.2698969179570267E-002, 8.2823552555632743E-002, 7.3139765360780429E-002,
      6.3819305996729803E-002, 5.5000931756790189E-002, 4.6767242937952071E-002,
      3.9166540734920929E-002, 3.2230799516104032E-002, 2.5975728678253356E-002,
      2.0398539039679763E-002, 1.5481539214729679E-002, 1.1197252690703319E-002,
      7.5112309421968976E-003, 4.3834710955235959E-003, 1.7700531216026730E-003,
      -3.7501733915976034E-004, -2.0983385200421876E-003,
      -3.4460892573046401E-003, -4.4632234248221754E-003,
      -5.1927820420139548E-003, -5.6753591527284510E-003,
      -5.9487386229897200E-003, -6.0476730056964996E-003,
      -6.0037714249869525E-003, -5.8454610152542776E-003,
      -5.5979908462344330E-003, -5.2835069506470119E-003,
      -4.9212628611336028E-003, -4.5278784030896796E-003,
      -4.1173854825626984E-003, -3.7011285585554169E-003,
      -3.2882268113373146E-003, -2.8867456953911310E-003,
      -2.5037338534440062E-003, -2.1427901412592624E-003,
      -1.8034258430109087E-003, -1.4879239647352142E-003,
      -1.2070965440010295E-003, -9.6534337899673093E-004,
      -7.3828062143524778E-004, -5.0212923771870239E-004,
      -3.0884957682706109E-004, -2.4525655082146565E-004,
      -2.1141541893125084E-004, 8.0118879554537623E-005, 5.1002343310981740E-004,
      2.5298926799180010E-004, -7.7799875809259262E-004,
      -4.0757125398815312E-004, 2.6895817890755605E-003, 3.3454344736318512E-003,
      -4.5730299065465880E-003, -1.0264948272560850E-002,
      8.1198879845256391E-003, 3.2531918686394358E-002, -4.6528912316376944E-003,
      -8.8810561877557856E-002, -3.1250113656430219E-002,
      2.2308979719173758E-001, 1.9561190800160422E-001, -5.0051653523851081E-001,
      -7.7488559344793084E-001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
