/*
 * helikopterExercise4a.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterExercise4a.mdl".
 *
 * Model version              : 1.64
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Tue Apr 28 12:54:45 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopterExercise4a.h"
#include "helikopterExercise4a_private.h"
#include <stdio.h>
#include "helikopterExercise4a_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopterExercise4a helikopterExercise4a_B;

/* Continuous states */
ContinuousStates_helikopterExer helikopterExercise4a_X;

/* Block states (auto storage) */
D_Work_helikopterExercise4a helikopterExercise4a_DWork;

/* Real-time model */
RT_MODEL_helikopterExercise4a helikopterExercise4a_M_;
RT_MODEL_helikopterExercise4a *helikopterExercise4a_M = &helikopterExercise4a_M_;

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
  helikopterExercise4a_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopterExercise4a_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Gain2;
  real_T rtb_Saturation_m;
  real_T rtb_Gain1;
  real_T rtb_Gain[6];
  real_T rtb_K_ed;
  real_T tmp[6];
  int32_T tmp_0;
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
    /* set solver stop time */
    if (!(helikopterExercise4a_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopterExercise4a_M->solverInfo,
                            ((helikopterExercise4a_M->Timing.clockTickH0 + 1) *
        helikopterExercise4a_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopterExercise4a_M->solverInfo,
                            ((helikopterExercise4a_M->Timing.clockTick0 + 1) *
        helikopterExercise4a_M->Timing.stepSize0 +
        helikopterExercise4a_M->Timing.clockTickH0 *
        helikopterExercise4a_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopterExercise4a_M)) {
    helikopterExercise4a_M->Timing.t[0] = rtsiGetT
      (&helikopterExercise4a_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise4a_B.VandringLavpass =
    helikopterExercise4a_P.VandringLavpass_C*
    helikopterExercise4a_X.VandringLavpass_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
      if (!(++helikopterExercise4a_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopterExercise4a_DWork.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise4a_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4a_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopterExercise4a_M->Timing.t[1];
          u[1] = helikopterExercise4a_B.VandringLavpass;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4a_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopterExercise4a_DWork.ToFile_IWORK.Count)*2)+1 >=
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

    /* S-Function Block: helikopterExercise4a/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder
        (helikopterExercise4a_DWork.HILInitialize_Card,
         helikopterExercise4a_P.HILReadEncoder_Channels, 3,
         &helikopterExercise4a_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 =
          helikopterExercise4a_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 =
          helikopterExercise4a_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 =
          helikopterExercise4a_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterExercise4a_B.KalibrerPitch =
      helikopterExercise4a_P.KalibrerPitch_Gain * rtb_HILReadEncoder_o2;

    /* ToFile: '<Root>/To File1' */
    if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
      if (!(++helikopterExercise4a_DWork.ToFile1_IWORK.Decimation % 1) &&
          (helikopterExercise4a_DWork.ToFile1_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise4a_DWork.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4a_DWork.ToFile1_IWORK.Decimation = 0;
          u[0] = helikopterExercise4a_M->Timing.t[1];
          u[1] = helikopterExercise4a_B.KalibrerPitch;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4a_M,
                              "Error writing to MAT-file pitch.mat");
            return;
          }

          if (((++helikopterExercise4a_DWork.ToFile1_IWORK.Count)*2)+1 >=
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
    helikopterExercise4a_B.KalibrerElev =
      helikopterExercise4a_P.KalibrerElev_Gain * rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopterExercise4a_B.Add = helikopterExercise4a_B.KalibrerElev +
      helikopterExercise4a_P.Constant_Value;

    /* ToFile: '<Root>/To File2' */
    if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
      if (!(++helikopterExercise4a_DWork.ToFile2_IWORK.Decimation % 1) &&
          (helikopterExercise4a_DWork.ToFile2_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise4a_DWork.ToFile2_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4a_DWork.ToFile2_IWORK.Decimation = 0;
          u[0] = helikopterExercise4a_M->Timing.t[1];
          u[1] = helikopterExercise4a_B.Add;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4a_M,
                              "Error writing to MAT-file elevation.mat");
            return;
          }

          if (((++helikopterExercise4a_DWork.ToFile2_IWORK.Count)*2)+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file elevation.mat.\n");
          }
        }
      }
    }
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopterExercise4a_X.Integrator_CSTATE >=
      helikopterExercise4a_P.Integrator_UpperSat ) {
    helikopterExercise4a_X.Integrator_CSTATE =
      helikopterExercise4a_P.Integrator_UpperSat;
  } else if (helikopterExercise4a_X.Integrator_CSTATE <=
             helikopterExercise4a_P.Integrator_LowerSat ) {
    helikopterExercise4a_X.Integrator_CSTATE =
      helikopterExercise4a_P.Integrator_LowerSat;
  }

  rtb_Saturation_m = helikopterExercise4a_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopterExercise4a_B.KalibrerVandring =
      helikopterExercise4a_P.KalibrerVandring_Gain * rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopterExercise4a_P.VandringDeriv_D*
    helikopterExercise4a_B.KalibrerVandring;
  rtb_VandringDeriv += helikopterExercise4a_P.VandringDeriv_C*
    helikopterExercise4a_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopterExercise4a_P.TransferFcn4_D*
    helikopterExercise4a_B.KalibrerPitch;
  rtb_Gain2 += helikopterExercise4a_P.TransferFcn4_C*
    helikopterExercise4a_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopterExercise4a_P.TransferFcn5_D*
    helikopterExercise4a_B.KalibrerElev;
  rtb_Gain1 += helikopterExercise4a_P.TransferFcn5_C*
    helikopterExercise4a_X.TransferFcn5_CSTATE;

  /* Gain: '<Root>/Gain' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtGainInport1'
   */
  tmp[0] = helikopterExercise4a_B.VandringLavpass;
  tmp[1] = rtb_VandringDeriv;
  tmp[2] = helikopterExercise4a_B.KalibrerPitch;
  tmp[3] = rtb_Gain2;
  tmp[4] = helikopterExercise4a_B.Add;
  tmp[5] = rtb_Gain1;
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_Gain[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain[tmp_0] += helikopterExercise4a_P.Gain_Gain[6 * tmp_1 + tmp_0] *
        tmp[tmp_1];
    }
  }

  /* Gain: '<S1>/K_ed' */
  rtb_K_ed = helikopterExercise4a_P.K_ed_Gain * rtb_Gain[5];

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise4a_DWork.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise4a_DWork.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex =
      helikopterExercise4a_DWork.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helikopterExercise4a_M->Timing.t[0];

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

    helikopterExercise4a_DWork.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

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
        pDataValues += 81;
      }
    }
  }

  /* Sum: '<S1>/Sum' */
  rtb_Gain1 -= rtb_Gain[4];

  /* Gain: '<S1>/K_ei' */
  helikopterExercise4a_B.K_ei = helikopterExercise4a_P.K_ei_Gain * rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_m = (helikopterExercise4a_P.K_ep_Gain * rtb_Gain1 +
                      rtb_Saturation_m) + rtb_K_ed;

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_m = rt_SATURATE(rtb_Saturation_m,
    helikopterExercise4a_P.Saturation_LowerSat,
    helikopterExercise4a_P.Saturation_UpperSat);

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise4a_DWork.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise4a_DWork.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex =
      helikopterExercise4a_DWork.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopterExercise4a_M->Timing.t[0];

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

    helikopterExercise4a_DWork.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

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
        pDataValues += 81;
      }
    }
  }

  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/K_pd'
   *  Gain: '<S3>/K_pp'
   *  Saturate: '<S3>/Saturation'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Gain1,
    helikopterExercise4a_P.Saturation_LowerSat_g,
    helikopterExercise4a_P.Saturation_UpperSat_i) - rtb_Gain[2]) *
    helikopterExercise4a_P.K_pp_Gain - helikopterExercise4a_P.K_pd_Gain *
    rtb_Gain[3];

  /* Gain: '<S4>/Gain2' incorporates:
   *  Sum: '<S4>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_m - rtb_Gain1) * helikopterExercise4a_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopterExercise4a_B.SatB = rt_SATURATE(rtb_Gain2,
    helikopterExercise4a_P.SatB_LowerSat, helikopterExercise4a_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
  }

  /* Gain: '<S4>/Gain1' incorporates:
   *  Sum: '<S4>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_m) * helikopterExercise4a_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopterExercise4a_B.Sat = rt_SATURATE(rtb_Gain1,
    helikopterExercise4a_P.Sat_LowerSat, helikopterExercise4a_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopterExercise4a/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopterExercise4a_DWork.HILWriteAnalog_Buffer[0] =
        helikopterExercise4a_B.SatB;
      helikopterExercise4a_DWork.HILWriteAnalog_Buffer[1] =
        helikopterExercise4a_B.Sat;
      result = hil_write_analog(helikopterExercise4a_DWork.HILInitialize_Card,
        helikopterExercise4a_P.HILWriteAnalog_Channels, 2,
        &helikopterExercise4a_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopterExercise4a_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
    rt_ertODEUpdateContinuousStates(&helikopterExercise4a_M->solverInfo);
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
  if (!(++helikopterExercise4a_M->Timing.clockTick0)) {
    ++helikopterExercise4a_M->Timing.clockTickH0;
  }

  helikopterExercise4a_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopterExercise4a_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopterExercise4a_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopterExercise4a_M->Timing.clockTick1)) {
      ++helikopterExercise4a_M->Timing.clockTickH1;
    }

    helikopterExercise4a_M->Timing.t[1] =
      helikopterExercise4a_M->Timing.clockTick1 *
      helikopterExercise4a_M->Timing.stepSize1 +
      helikopterExercise4a_M->Timing.clockTickH1 *
      helikopterExercise4a_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopterExercise4a_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->VandringLavpass_CSTATE =
      helikopterExercise4a_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->VandringLavpass_CSTATE +=
      (helikopterExercise4a_P.VandringLavpass_A)*
      helikopterExercise4a_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopterExercise4a_X.Integrator_CSTATE <=
            helikopterExercise4a_P.Integrator_LowerSat );
    usat = ( helikopterExercise4a_X.Integrator_CSTATE >=
            helikopterExercise4a_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopterExercise4a_B.K_ei > 0)) ||
        (usat && (helikopterExercise4a_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopterExer *)
        helikopterExercise4a_M->ModelData.derivs)->Integrator_CSTATE =
        helikopterExercise4a_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopterExer *)
        helikopterExercise4a_M->ModelData.derivs)->Integrator_CSTATE = 0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->VandringDeriv_CSTATE =
      helikopterExercise4a_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->VandringDeriv_CSTATE +=
      (helikopterExercise4a_P.VandringDeriv_A)*
      helikopterExercise4a_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->TransferFcn4_CSTATE =
      helikopterExercise4a_B.KalibrerPitch;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->TransferFcn4_CSTATE +=
      (helikopterExercise4a_P.TransferFcn4_A)*
      helikopterExercise4a_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->TransferFcn5_CSTATE =
      helikopterExercise4a_B.KalibrerElev;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4a_M->ModelData.derivs)->TransferFcn5_CSTATE +=
      (helikopterExercise4a_P.TransferFcn5_A)*
      helikopterExercise4a_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopterExercise4a_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopterExercise4a_P.Integrator_UpperSat = rtInf;
  helikopterExercise4a_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopterExercise4a_M, 0,
                sizeof(RT_MODEL_helikopterExercise4a));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopterExercise4a_M->solverInfo,
                          &helikopterExercise4a_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopterExercise4a_M->solverInfo, &rtmGetTPtr
                (helikopterExercise4a_M));
    rtsiSetStepSizePtr(&helikopterExercise4a_M->solverInfo,
                       &helikopterExercise4a_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopterExercise4a_M->solverInfo,
                 &helikopterExercise4a_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopterExercise4a_M->solverInfo,
                         &helikopterExercise4a_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopterExercise4a_M->solverInfo,
      &helikopterExercise4a_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopterExercise4a_M->solverInfo,
                          (&rtmGetErrorStatus(helikopterExercise4a_M)));
    rtsiSetRTModelPtr(&helikopterExercise4a_M->solverInfo,
                      helikopterExercise4a_M);
  }

  rtsiSetSimTimeStep(&helikopterExercise4a_M->solverInfo, MAJOR_TIME_STEP);
  helikopterExercise4a_M->ModelData.intgData.f[0] =
    helikopterExercise4a_M->ModelData.odeF[0];
  helikopterExercise4a_M->ModelData.contStates = ((real_T *)
    &helikopterExercise4a_X);
  rtsiSetSolverData(&helikopterExercise4a_M->solverInfo, (void *)
                    &helikopterExercise4a_M->ModelData.intgData);
  rtsiSetSolverName(&helikopterExercise4a_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopterExercise4a_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopterExercise4a_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopterExercise4a_M->Timing.sampleTimes =
      (&helikopterExercise4a_M->Timing.sampleTimesArray[0]);
    helikopterExercise4a_M->Timing.offsetTimes =
      (&helikopterExercise4a_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopterExercise4a_M->Timing.sampleTimes[0] = (0.0);
    helikopterExercise4a_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopterExercise4a_M->Timing.offsetTimes[0] = (0.0);
    helikopterExercise4a_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopterExercise4a_M, &helikopterExercise4a_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopterExercise4a_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopterExercise4a_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopterExercise4a_M, -1);
  helikopterExercise4a_M->Timing.stepSize0 = 0.001;
  helikopterExercise4a_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopterExercise4a_M->Sizes.checksums[0] = (4145194891U);
  helikopterExercise4a_M->Sizes.checksums[1] = (3570954371U);
  helikopterExercise4a_M->Sizes.checksums[2] = (4202034545U);
  helikopterExercise4a_M->Sizes.checksums[3] = (2600162562U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopterExercise4a_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopterExercise4a_M->extModeInfo,
      &helikopterExercise4a_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopterExercise4a_M->extModeInfo,
                        helikopterExercise4a_M->Sizes.checksums);
    rteiSetTPtr(helikopterExercise4a_M->extModeInfo, rtmGetTPtr
                (helikopterExercise4a_M));
  }

  helikopterExercise4a_M->solverInfoPtr = (&helikopterExercise4a_M->solverInfo);
  helikopterExercise4a_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopterExercise4a_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopterExercise4a_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopterExercise4a_M->ModelData.blockIO = ((void *) &helikopterExercise4a_B);

  {
    helikopterExercise4a_B.VandringLavpass = 0.0;
    helikopterExercise4a_B.KalibrerPitch = 0.0;
    helikopterExercise4a_B.KalibrerElev = 0.0;
    helikopterExercise4a_B.Add = 0.0;
    helikopterExercise4a_B.KalibrerVandring = 0.0;
    helikopterExercise4a_B.K_ei = 0.0;
    helikopterExercise4a_B.SatB = 0.0;
    helikopterExercise4a_B.Sat = 0.0;
  }

  /* parameters */
  helikopterExercise4a_M->ModelData.defaultParam = ((real_T *)
    &helikopterExercise4a_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopterExercise4a_X;
    helikopterExercise4a_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopterExercise4a_X, 0,
                  sizeof(ContinuousStates_helikopterExer));
  }

  /* states (dwork) */
  helikopterExercise4a_M->Work.dwork = ((void *) &helikopterExercise4a_DWork);
  (void) memset((void *)&helikopterExercise4a_DWork, 0,
                sizeof(D_Work_helikopterExercise4a));

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise4a_DWork.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise4a_DWork.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  helikopterExercise4a_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopterExercise4a_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopterExercise4a_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopterExercise4a_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopterExercise4a_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopterExercise4a_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise4a/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopterExercise4a_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopterExercise4a_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopterExercise4a_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopterExercise4a_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopterExercise4a_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopterExercise4a_P.HILInitialize_AOExit && is_switching)) {
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise4a_P.HILInitialize_AOFinal;
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise4a_P.HILInitialize_AOFinal;
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise4a_P.HILInitialize_AOFinal;
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise4a_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopterExercise4a_DWork.HILInitialize_Card,
        &helikopterExercise4a_P.HILInitialize_AOChannels[0], 4U,
        &helikopterExercise4a_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
      }
    }

    hil_close(helikopterExercise4a_DWork.HILInitialize_Card);
    helikopterExercise4a_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopterExercise4a_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4a_DWork.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error writing header for ans to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      helikopterExercise4a_DWork.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) helikopterExercise4a_DWork.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "pitch.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error reopening MAT-file pitch.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4a_DWork.ToFile1_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error writing header for ans to MAT-file pitch.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      helikopterExercise4a_DWork.ToFile1_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File2' */
  {
    FILE *fp = (FILE *) helikopterExercise4a_DWork.ToFile2_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "elevation.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error closing MAT-file elevation.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error reopening MAT-file elevation.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4a_DWork.ToFile2_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error writing header for ans to MAT-file elevation.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4a_M,
                          "Error closing MAT-file elevation.mat");
        return;
      }

      helikopterExercise4a_DWork.ToFile2_PWORK.FilePtr = (NULL);
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
  helikopterExercise4a_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopterExercise4a_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopterExercise4a_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopterExercise4a_M->Sizes.numY = (0);/* Number of model outputs */
  helikopterExercise4a_M->Sizes.numU = (0);/* Number of model inputs */
  helikopterExercise4a_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopterExercise4a_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopterExercise4a_M->Sizes.numBlocks = (46);/* Number of blocks */
  helikopterExercise4a_M->Sizes.numBlockIO = (8);/* Number of block outputs */
  helikopterExercise4a_M->Sizes.numBlockPrms = (148);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise4a_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopterExercise4a_X.Integrator_CSTATE =
    helikopterExercise4a_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopterExercise4a_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopterExercise4a_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopterExercise4a_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise4a/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("sensoray_model_626", "0",
                      &helikopterExercise4a_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopterExercise4a_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopterExercise4a_P.HILInitialize_CKPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_ClockModes =
          &helikopterExercise4a_DWork.HILInitialize_ClockModes[0];
        for (i1=0; i1 < 6; i1++) {
          dw_ClockModes[i1] = helikopterExercise4a_P.HILInitialize_CKModes;
        }
      }

      result = hil_set_clock_mode(helikopterExercise4a_DWork.HILInitialize_Card,
                                  (t_clock *)
        &helikopterExercise4a_P.HILInitialize_CKChannels[0], 6U, (t_clock_mode *)
        &helikopterExercise4a_DWork.HILInitialize_ClockModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopterExercise4a_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
      return;
    }

    if ((helikopterExercise4a_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopterExercise4a_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helikopterExercise4a_DWork.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMinimums[i1] = helikopterExercise4a_P.HILInitialize_AILow;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helikopterExercise4a_DWork.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMaximums[i1] = helikopterExercise4a_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helikopterExercise4a_DWork.HILInitialize_Card,
         &helikopterExercise4a_P.HILInitialize_AIChannels[0], 16U,
         &helikopterExercise4a_DWork.HILInitialize_AIMinimums[0],
         &helikopterExercise4a_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise4a_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopterExercise4a_P.HILInitialize_AOPEnter && is_switching)) {
      helikopterExercise4a_DWork.HILInitialize_AOMinimums[0] =
        helikopterExercise4a_P.HILInitialize_AOLow;
      helikopterExercise4a_DWork.HILInitialize_AOMinimums[1] =
        helikopterExercise4a_P.HILInitialize_AOLow;
      helikopterExercise4a_DWork.HILInitialize_AOMinimums[2] =
        helikopterExercise4a_P.HILInitialize_AOLow;
      helikopterExercise4a_DWork.HILInitialize_AOMinimums[3] =
        helikopterExercise4a_P.HILInitialize_AOLow;
      helikopterExercise4a_DWork.HILInitialize_AOMaximums[0] =
        helikopterExercise4a_P.HILInitialize_AOHigh;
      helikopterExercise4a_DWork.HILInitialize_AOMaximums[1] =
        helikopterExercise4a_P.HILInitialize_AOHigh;
      helikopterExercise4a_DWork.HILInitialize_AOMaximums[2] =
        helikopterExercise4a_P.HILInitialize_AOHigh;
      helikopterExercise4a_DWork.HILInitialize_AOMaximums[3] =
        helikopterExercise4a_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopterExercise4a_DWork.HILInitialize_Card,
         &helikopterExercise4a_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4a_DWork.HILInitialize_AOMinimums[0],
         &helikopterExercise4a_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise4a_P.HILInitialize_AOStart && !is_switching) ||
        (helikopterExercise4a_P.HILInitialize_AOEnter && is_switching)) {
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise4a_P.HILInitialize_AOInitial;
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise4a_P.HILInitialize_AOInitial;
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise4a_P.HILInitialize_AOInitial;
      helikopterExercise4a_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise4a_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopterExercise4a_DWork.HILInitialize_Card,
        &helikopterExercise4a_P.HILInitialize_AOChannels[0], 4U,
        &helikopterExercise4a_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4a_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4a_M,
                        "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4a_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopterExercise4a_DWork.ToFile_IWORK.Count = 0;
    helikopterExercise4a_DWork.ToFile_IWORK.Decimation = -1;
    helikopterExercise4a_DWork.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File1' */
  {
    const char *fileName = "pitch.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4a_M,
                        "Error creating .mat file pitch.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4a_M,
                        "Error writing mat file header to file pitch.mat");
      return;
    }

    helikopterExercise4a_DWork.ToFile1_IWORK.Count = 0;
    helikopterExercise4a_DWork.ToFile1_IWORK.Decimation = -1;
    helikopterExercise4a_DWork.ToFile1_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File2' */
  {
    const char *fileName = "elevation.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4a_M,
                        "Error creating .mat file elevation.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4a_M,
                        "Error writing mat file header to file elevation.mat");
      return;
    }

    helikopterExercise4a_DWork.ToFile2_IWORK.Count = 0;
    helikopterExercise4a_DWork.ToFile2_IWORK.Decimation = -1;
    helikopterExercise4a_DWork.ToFile2_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
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
      8.5719423857623291E-003, 1.2221744528993380E-002, 1.6638702459657106E-002,
      2.3031297902052886E-002, 3.1832678266465832E-002, 4.3668613460816942E-002,
      5.9537225181344577E-002, 8.0752100009053890E-002, 1.0825027178731837E-001,
      1.4304083525633746E-001, 1.8475954244149273E-001, 2.3060586368680336E-001,
      2.7151577427895607E-001, 2.8607575066473145E-001, 2.2748545611083540E-001,
      -2.6620032527911286E-005, -6.0600428787701833E-005,
      -4.0772964153704746E-007, 2.7249278319789954E-005,
      -4.2782725177049062E-007, 2.7888082637255269E-005,
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

    helikopterExercise4a_DWork.FromWorkspace1_PWORK.TimePtr = (void *)
      pTimeValues;
    helikopterExercise4a_DWork.FromWorkspace1_PWORK.DataPtr = (void *)
      pDataValues;
    helikopterExercise4a_DWork.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
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
      0.0 } ;

    helikopterExercise4a_DWork.FromWorkspace_PWORK.TimePtr = (void *)
      pTimeValues;
    helikopterExercise4a_DWork.FromWorkspace_PWORK.DataPtr = (void *)
      pDataValues;
    helikopterExercise4a_DWork.FromWorkspace_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopterExercise4a_terminate();
}

RT_MODEL_helikopterExercise4a *helikopterExercise4a(void)
{
  helikopterExercise4a_initialize(1);
  return helikopterExercise4a_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
