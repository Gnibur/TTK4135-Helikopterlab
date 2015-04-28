/*
 * helikopterExercise2.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterExercise2.mdl".
 *
 * Model version              : 1.63
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Tue Apr 28 12:16:02 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopterExercise2.h"
#include "helikopterExercise2_private.h"
#include <stdio.h>
#include "helikopterExercise2_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopterExercise2 helikopterExercise2_B;

/* Continuous states */
ContinuousStates_helikopterExer helikopterExercise2_X;

/* Block states (auto storage) */
D_Work_helikopterExercise2 helikopterExercise2_DWork;

/* Real-time model */
RT_MODEL_helikopterExercise2 helikopterExercise2_M_;
RT_MODEL_helikopterExercise2 *helikopterExercise2_M = &helikopterExercise2_M_;

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
  helikopterExercise2_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopterExercise2_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Gain2;
  real_T rtb_Saturation_i;
  real_T rtb_Gain1;
  real_T rtb_Gain[6];
  real_T tmp[6];
  int32_T tmp_0;
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
    /* set solver stop time */
    if (!(helikopterExercise2_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopterExercise2_M->solverInfo,
                            ((helikopterExercise2_M->Timing.clockTickH0 + 1) *
        helikopterExercise2_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopterExercise2_M->solverInfo,
                            ((helikopterExercise2_M->Timing.clockTick0 + 1) *
        helikopterExercise2_M->Timing.stepSize0 +
        helikopterExercise2_M->Timing.clockTickH0 *
        helikopterExercise2_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopterExercise2_M)) {
    helikopterExercise2_M->Timing.t[0] = rtsiGetT
      (&helikopterExercise2_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise2_B.VandringLavpass =
    helikopterExercise2_P.VandringLavpass_C*
    helikopterExercise2_X.VandringLavpass_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
      if (!(++helikopterExercise2_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopterExercise2_DWork.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise2_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise2_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopterExercise2_M->Timing.t[1];
          u[1] = helikopterExercise2_B.VandringLavpass;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise2_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopterExercise2_DWork.ToFile_IWORK.Count)*2)+1 >=
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

    /* S-Function Block: helikopterExercise2/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder
        (helikopterExercise2_DWork.HILInitialize_Card,
         helikopterExercise2_P.HILReadEncoder_Channels, 3,
         &helikopterExercise2_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 = helikopterExercise2_DWork.HILReadEncoder_Buffer
          [0];
        rtb_HILReadEncoder_o2 = helikopterExercise2_DWork.HILReadEncoder_Buffer
          [1];
        rtb_HILReadEncoder_o3 = helikopterExercise2_DWork.HILReadEncoder_Buffer
          [2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterExercise2_B.KalibrerPitch =
      helikopterExercise2_P.KalibrerPitch_Gain * rtb_HILReadEncoder_o2;

    /* ToFile: '<Root>/To File1' */
    if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
      if (!(++helikopterExercise2_DWork.ToFile1_IWORK.Decimation % 1) &&
          (helikopterExercise2_DWork.ToFile1_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) helikopterExercise2_DWork.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise2_DWork.ToFile1_IWORK.Decimation = 0;
          u[0] = helikopterExercise2_M->Timing.t[1];
          u[1] = helikopterExercise2_B.KalibrerPitch;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise2_M,
                              "Error writing to MAT-file pitch.mat");
            return;
          }

          if (((++helikopterExercise2_DWork.ToFile1_IWORK.Count)*2)+1 >=
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
    helikopterExercise2_B.KalibrerElev = helikopterExercise2_P.KalibrerElev_Gain
      * rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopterExercise2_B.Add = helikopterExercise2_B.KalibrerElev +
      helikopterExercise2_P.Constant_Value;
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopterExercise2_X.Integrator_CSTATE >=
      helikopterExercise2_P.Integrator_UpperSat ) {
    helikopterExercise2_X.Integrator_CSTATE =
      helikopterExercise2_P.Integrator_UpperSat;
  } else if (helikopterExercise2_X.Integrator_CSTATE <=
             helikopterExercise2_P.Integrator_LowerSat ) {
    helikopterExercise2_X.Integrator_CSTATE =
      helikopterExercise2_P.Integrator_LowerSat;
  }

  rtb_Saturation_i = helikopterExercise2_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopterExercise2_B.KalibrerVandring =
      helikopterExercise2_P.KalibrerVandring_Gain * rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopterExercise2_P.VandringDeriv_D*
    helikopterExercise2_B.KalibrerVandring;
  rtb_VandringDeriv += helikopterExercise2_P.VandringDeriv_C*
    helikopterExercise2_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopterExercise2_P.TransferFcn4_D*
    helikopterExercise2_B.KalibrerPitch;
  rtb_Gain2 += helikopterExercise2_P.TransferFcn4_C*
    helikopterExercise2_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopterExercise2_P.TransferFcn5_D*
    helikopterExercise2_B.KalibrerElev;
  rtb_Gain1 += helikopterExercise2_P.TransferFcn5_C*
    helikopterExercise2_X.TransferFcn5_CSTATE;

  /* Gain: '<Root>/Gain' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtGainInport1'
   */
  tmp[0] = helikopterExercise2_B.VandringLavpass;
  tmp[1] = rtb_VandringDeriv;
  tmp[2] = helikopterExercise2_B.KalibrerPitch;
  tmp[3] = rtb_Gain2;
  tmp[4] = helikopterExercise2_B.Add;
  tmp[5] = rtb_Gain1;
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_Gain[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain[tmp_0] += helikopterExercise2_P.Gain_Gain[6 * tmp_1 + tmp_0] *
        tmp[tmp_1];
    }
  }

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<Root>/elevation'
   */
  rtb_Gain1 = helikopterExercise2_P.elevation_Value - rtb_Gain[4];

  /* Gain: '<S1>/K_ei' */
  helikopterExercise2_B.K_ei = helikopterExercise2_P.K_ei_Gain * rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ed'
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_i = (helikopterExercise2_P.K_ep_Gain * rtb_Gain1 +
                      rtb_Saturation_i) + helikopterExercise2_P.K_ed_Gain *
    rtb_Gain[5];

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_i = rt_SATURATE(rtb_Saturation_i,
    helikopterExercise2_P.Saturation_LowerSat,
    helikopterExercise2_P.Saturation_UpperSat);

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise2_DWork.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise2_DWork.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex =
      helikopterExercise2_DWork.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopterExercise2_M->Timing.t[0];

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

    helikopterExercise2_DWork.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

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

  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/K_pd'
   *  Gain: '<S3>/K_pp'
   *  Saturate: '<S3>/Saturation'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Gain1,
    helikopterExercise2_P.Saturation_LowerSat_l,
    helikopterExercise2_P.Saturation_UpperSat_h) - rtb_Gain[2]) *
    helikopterExercise2_P.K_pp_Gain - helikopterExercise2_P.K_pd_Gain *
    rtb_Gain[3];

  /* Gain: '<S4>/Gain2' incorporates:
   *  Sum: '<S4>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_i - rtb_Gain1) * helikopterExercise2_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopterExercise2_B.SatB = rt_SATURATE(rtb_Gain2,
    helikopterExercise2_P.SatB_LowerSat, helikopterExercise2_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
  }

  /* Gain: '<S4>/Gain1' incorporates:
   *  Sum: '<S4>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_i) * helikopterExercise2_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopterExercise2_B.Sat = rt_SATURATE(rtb_Gain1,
    helikopterExercise2_P.Sat_LowerSat, helikopterExercise2_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopterExercise2/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopterExercise2_DWork.HILWriteAnalog_Buffer[0] =
        helikopterExercise2_B.SatB;
      helikopterExercise2_DWork.HILWriteAnalog_Buffer[1] =
        helikopterExercise2_B.Sat;
      result = hil_write_analog(helikopterExercise2_DWork.HILInitialize_Card,
        helikopterExercise2_P.HILWriteAnalog_Channels, 2,
        &helikopterExercise2_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopterExercise2_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
    rt_ertODEUpdateContinuousStates(&helikopterExercise2_M->solverInfo);
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
  if (!(++helikopterExercise2_M->Timing.clockTick0)) {
    ++helikopterExercise2_M->Timing.clockTickH0;
  }

  helikopterExercise2_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopterExercise2_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopterExercise2_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopterExercise2_M->Timing.clockTick1)) {
      ++helikopterExercise2_M->Timing.clockTickH1;
    }

    helikopterExercise2_M->Timing.t[1] =
      helikopterExercise2_M->Timing.clockTick1 *
      helikopterExercise2_M->Timing.stepSize1 +
      helikopterExercise2_M->Timing.clockTickH1 *
      helikopterExercise2_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopterExercise2_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->VandringLavpass_CSTATE = helikopterExercise2_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->VandringLavpass_CSTATE += (helikopterExercise2_P.VandringLavpass_A)*
      helikopterExercise2_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopterExercise2_X.Integrator_CSTATE <=
            helikopterExercise2_P.Integrator_LowerSat );
    usat = ( helikopterExercise2_X.Integrator_CSTATE >=
            helikopterExercise2_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopterExercise2_B.K_ei > 0)) ||
        (usat && (helikopterExercise2_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopterExer *)
        helikopterExercise2_M->ModelData.derivs)->Integrator_CSTATE =
        helikopterExercise2_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopterExer *)
        helikopterExercise2_M->ModelData.derivs)->Integrator_CSTATE = 0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->VandringDeriv_CSTATE = helikopterExercise2_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->VandringDeriv_CSTATE += (helikopterExercise2_P.VandringDeriv_A)*
      helikopterExercise2_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->TransferFcn4_CSTATE = helikopterExercise2_B.KalibrerPitch;
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->TransferFcn4_CSTATE += (helikopterExercise2_P.TransferFcn4_A)*
      helikopterExercise2_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->TransferFcn5_CSTATE = helikopterExercise2_B.KalibrerElev;
    ((StateDerivatives_helikopterExer *) helikopterExercise2_M->ModelData.derivs)
      ->TransferFcn5_CSTATE += (helikopterExercise2_P.TransferFcn5_A)*
      helikopterExercise2_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopterExercise2_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopterExercise2_P.Integrator_UpperSat = rtInf;
  helikopterExercise2_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopterExercise2_M, 0,
                sizeof(RT_MODEL_helikopterExercise2));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopterExercise2_M->solverInfo,
                          &helikopterExercise2_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopterExercise2_M->solverInfo, &rtmGetTPtr
                (helikopterExercise2_M));
    rtsiSetStepSizePtr(&helikopterExercise2_M->solverInfo,
                       &helikopterExercise2_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopterExercise2_M->solverInfo,
                 &helikopterExercise2_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopterExercise2_M->solverInfo,
                         &helikopterExercise2_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopterExercise2_M->solverInfo,
      &helikopterExercise2_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopterExercise2_M->solverInfo,
                          (&rtmGetErrorStatus(helikopterExercise2_M)));
    rtsiSetRTModelPtr(&helikopterExercise2_M->solverInfo, helikopterExercise2_M);
  }

  rtsiSetSimTimeStep(&helikopterExercise2_M->solverInfo, MAJOR_TIME_STEP);
  helikopterExercise2_M->ModelData.intgData.f[0] =
    helikopterExercise2_M->ModelData.odeF[0];
  helikopterExercise2_M->ModelData.contStates = ((real_T *)
    &helikopterExercise2_X);
  rtsiSetSolverData(&helikopterExercise2_M->solverInfo, (void *)
                    &helikopterExercise2_M->ModelData.intgData);
  rtsiSetSolverName(&helikopterExercise2_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopterExercise2_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopterExercise2_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopterExercise2_M->Timing.sampleTimes =
      (&helikopterExercise2_M->Timing.sampleTimesArray[0]);
    helikopterExercise2_M->Timing.offsetTimes =
      (&helikopterExercise2_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopterExercise2_M->Timing.sampleTimes[0] = (0.0);
    helikopterExercise2_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopterExercise2_M->Timing.offsetTimes[0] = (0.0);
    helikopterExercise2_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopterExercise2_M, &helikopterExercise2_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopterExercise2_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopterExercise2_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopterExercise2_M, -1);
  helikopterExercise2_M->Timing.stepSize0 = 0.001;
  helikopterExercise2_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopterExercise2_M->Sizes.checksums[0] = (1888663190U);
  helikopterExercise2_M->Sizes.checksums[1] = (3805369506U);
  helikopterExercise2_M->Sizes.checksums[2] = (3515205500U);
  helikopterExercise2_M->Sizes.checksums[3] = (1706016471U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopterExercise2_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopterExercise2_M->extModeInfo,
      &helikopterExercise2_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopterExercise2_M->extModeInfo,
                        helikopterExercise2_M->Sizes.checksums);
    rteiSetTPtr(helikopterExercise2_M->extModeInfo, rtmGetTPtr
                (helikopterExercise2_M));
  }

  helikopterExercise2_M->solverInfoPtr = (&helikopterExercise2_M->solverInfo);
  helikopterExercise2_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopterExercise2_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopterExercise2_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopterExercise2_M->ModelData.blockIO = ((void *) &helikopterExercise2_B);

  {
    helikopterExercise2_B.VandringLavpass = 0.0;
    helikopterExercise2_B.KalibrerPitch = 0.0;
    helikopterExercise2_B.KalibrerElev = 0.0;
    helikopterExercise2_B.Add = 0.0;
    helikopterExercise2_B.KalibrerVandring = 0.0;
    helikopterExercise2_B.K_ei = 0.0;
    helikopterExercise2_B.SatB = 0.0;
    helikopterExercise2_B.Sat = 0.0;
  }

  /* parameters */
  helikopterExercise2_M->ModelData.defaultParam = ((real_T *)
    &helikopterExercise2_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopterExercise2_X;
    helikopterExercise2_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopterExercise2_X, 0,
                  sizeof(ContinuousStates_helikopterExer));
  }

  /* states (dwork) */
  helikopterExercise2_M->Work.dwork = ((void *) &helikopterExercise2_DWork);
  (void) memset((void *)&helikopterExercise2_DWork, 0,
                sizeof(D_Work_helikopterExercise2));

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise2_DWork.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise2_DWork.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  helikopterExercise2_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopterExercise2_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopterExercise2_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopterExercise2_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopterExercise2_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopterExercise2_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise2/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopterExercise2_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopterExercise2_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopterExercise2_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopterExercise2_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopterExercise2_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopterExercise2_P.HILInitialize_AOExit && is_switching)) {
      helikopterExercise2_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise2_P.HILInitialize_AOFinal;
      helikopterExercise2_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise2_P.HILInitialize_AOFinal;
      helikopterExercise2_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise2_P.HILInitialize_AOFinal;
      helikopterExercise2_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise2_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopterExercise2_DWork.HILInitialize_Card,
        &helikopterExercise2_P.HILInitialize_AOChannels[0], 4U,
        &helikopterExercise2_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
      }
    }

    hil_close(helikopterExercise2_DWork.HILInitialize_Card);
    helikopterExercise2_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopterExercise2_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise2_DWork.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error writing header for ans to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      helikopterExercise2_DWork.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) helikopterExercise2_DWork.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "pitch.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error reopening MAT-file pitch.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise2_DWork.ToFile1_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error writing header for ans to MAT-file pitch.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise2_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      helikopterExercise2_DWork.ToFile1_PWORK.FilePtr = (NULL);
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
  helikopterExercise2_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopterExercise2_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopterExercise2_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopterExercise2_M->Sizes.numY = (0);/* Number of model outputs */
  helikopterExercise2_M->Sizes.numU = (0);/* Number of model inputs */
  helikopterExercise2_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopterExercise2_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopterExercise2_M->Sizes.numBlocks = (45);/* Number of blocks */
  helikopterExercise2_M->Sizes.numBlockIO = (8);/* Number of block outputs */
  helikopterExercise2_M->Sizes.numBlockPrms = (149);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise2_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopterExercise2_X.Integrator_CSTATE = helikopterExercise2_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopterExercise2_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopterExercise2_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopterExercise2_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise2/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("sensoray_model_626", "0",
                      &helikopterExercise2_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopterExercise2_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopterExercise2_P.HILInitialize_CKPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_ClockModes =
          &helikopterExercise2_DWork.HILInitialize_ClockModes[0];
        for (i1=0; i1 < 6; i1++) {
          dw_ClockModes[i1] = helikopterExercise2_P.HILInitialize_CKModes;
        }
      }

      result = hil_set_clock_mode(helikopterExercise2_DWork.HILInitialize_Card,
                                  (t_clock *)
        &helikopterExercise2_P.HILInitialize_CKChannels[0], 6U, (t_clock_mode *)
        &helikopterExercise2_DWork.HILInitialize_ClockModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopterExercise2_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
      return;
    }

    if ((helikopterExercise2_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopterExercise2_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helikopterExercise2_DWork.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMinimums[i1] = helikopterExercise2_P.HILInitialize_AILow;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helikopterExercise2_DWork.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMaximums[i1] = helikopterExercise2_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helikopterExercise2_DWork.HILInitialize_Card,
         &helikopterExercise2_P.HILInitialize_AIChannels[0], 16U,
         &helikopterExercise2_DWork.HILInitialize_AIMinimums[0],
         &helikopterExercise2_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise2_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopterExercise2_P.HILInitialize_AOPEnter && is_switching)) {
      helikopterExercise2_DWork.HILInitialize_AOMinimums[0] =
        helikopterExercise2_P.HILInitialize_AOLow;
      helikopterExercise2_DWork.HILInitialize_AOMinimums[1] =
        helikopterExercise2_P.HILInitialize_AOLow;
      helikopterExercise2_DWork.HILInitialize_AOMinimums[2] =
        helikopterExercise2_P.HILInitialize_AOLow;
      helikopterExercise2_DWork.HILInitialize_AOMinimums[3] =
        helikopterExercise2_P.HILInitialize_AOLow;
      helikopterExercise2_DWork.HILInitialize_AOMaximums[0] =
        helikopterExercise2_P.HILInitialize_AOHigh;
      helikopterExercise2_DWork.HILInitialize_AOMaximums[1] =
        helikopterExercise2_P.HILInitialize_AOHigh;
      helikopterExercise2_DWork.HILInitialize_AOMaximums[2] =
        helikopterExercise2_P.HILInitialize_AOHigh;
      helikopterExercise2_DWork.HILInitialize_AOMaximums[3] =
        helikopterExercise2_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopterExercise2_DWork.HILInitialize_Card,
         &helikopterExercise2_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise2_DWork.HILInitialize_AOMinimums[0],
         &helikopterExercise2_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise2_P.HILInitialize_AOStart && !is_switching) ||
        (helikopterExercise2_P.HILInitialize_AOEnter && is_switching)) {
      helikopterExercise2_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise2_P.HILInitialize_AOInitial;
      helikopterExercise2_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise2_P.HILInitialize_AOInitial;
      helikopterExercise2_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise2_P.HILInitialize_AOInitial;
      helikopterExercise2_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise2_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopterExercise2_DWork.HILInitialize_Card,
        &helikopterExercise2_P.HILInitialize_AOChannels[0], 4U,
        &helikopterExercise2_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise2_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise2_M,
                        "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise2_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopterExercise2_DWork.ToFile_IWORK.Count = 0;
    helikopterExercise2_DWork.ToFile_IWORK.Decimation = -1;
    helikopterExercise2_DWork.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File1' */
  {
    const char *fileName = "pitch.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise2_M,
                        "Error creating .mat file pitch.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise2_M,
                        "Error writing mat file header to file pitch.mat");
      return;
    }

    helikopterExercise2_DWork.ToFile1_IWORK.Count = 0;
    helikopterExercise2_DWork.ToFile1_IWORK.Decimation = -1;
    helikopterExercise2_DWork.ToFile1_PWORK.FilePtr = fp;
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

    helikopterExercise2_DWork.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise2_DWork.FromWorkspace_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise2_DWork.FromWorkspace_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopterExercise2_terminate();
}

RT_MODEL_helikopterExercise2 *helikopterExercise2(void)
{
  helikopterExercise2_initialize(1);
  return helikopterExercise2_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
