/*
 * helikopterExercise4_feedback.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterExercise4_feedback.mdl".
 *
 * Model version              : 1.75
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Tue Apr 28 13:08:07 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopterExercise4_feedback.h"
#include "helikopterExercise4_feedback_private.h"
#include <stdio.h>
#include "helikopterExercise4_feedback_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopterExercise4_fee helikopterExercise4_feedback_B;

/* Continuous states */
ContinuousStates_helikopterExer helikopterExercise4_feedback_X;

/* Block states (auto storage) */
D_Work_helikopterExercise4_feed helikopterExercise4_feedb_DWork;

/* Real-time model */
RT_MODEL_helikopterExercise4_fe helikopterExercise4_feedback_M_;
RT_MODEL_helikopterExercise4_fe *helikopterExercise4_feedback_M =
  &helikopterExercise4_feedback_M_;

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
  helikopterExercise4_feedback_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopterExercise4_feedback_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Gain2;
  real_T rtb_Sum[2];
  real_T rtb_Saturation_g;
  real_T rtb_Sum1_p[6];
  real_T rtb_Gain1;
  real_T rtb_degrad[6];
  real_T rtb_K_ed;
  int32_T tmp;
  real_T tmp_0[2];
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
    /* set solver stop time */
    if (!(helikopterExercise4_feedback_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopterExercise4_feedback_M->solverInfo,
                            ((helikopterExercise4_feedback_M->Timing.clockTickH0
        + 1) * helikopterExercise4_feedback_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopterExercise4_feedback_M->solverInfo,
                            ((helikopterExercise4_feedback_M->Timing.clockTick0
        + 1) * helikopterExercise4_feedback_M->Timing.stepSize0 +
        helikopterExercise4_feedback_M->Timing.clockTickH0 *
        helikopterExercise4_feedback_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopterExercise4_feedback_M)) {
    helikopterExercise4_feedback_M->Timing.t[0] = rtsiGetT
      (&helikopterExercise4_feedback_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise4_feedback_B.VandringLavpass =
    helikopterExercise4_feedback_P.VandringLavpass_C*
    helikopterExercise4_feedback_X.VandringLavpass_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
    /* ToFile: '<Root>/To File' */
    if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
      if (!(++helikopterExercise4_feedb_DWork.ToFile_IWORK.Decimation % 1) &&
          (helikopterExercise4_feedb_DWork.ToFile_IWORK.Count*2)+1 < 100000000 )
      {
        FILE *fp = (FILE *) helikopterExercise4_feedb_DWork.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4_feedb_DWork.ToFile_IWORK.Decimation = 0;
          u[0] = helikopterExercise4_feedback_M->Timing.t[1];
          u[1] = helikopterExercise4_feedback_B.VandringLavpass;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4_feedback_M,
                              "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++helikopterExercise4_feedb_DWork.ToFile_IWORK.Count)*2)+1 >=
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

    /* S-Function Block: helikopterExercise4_feedback/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder
        (helikopterExercise4_feedb_DWork.HILInitialize_Card,
         helikopterExercise4_feedback_P.HILReadEncoder_Channels, 3,
         &helikopterExercise4_feedb_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 =
          helikopterExercise4_feedb_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 =
          helikopterExercise4_feedb_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 =
          helikopterExercise4_feedb_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopterExercise4_feedback_B.KalibrerElev =
      helikopterExercise4_feedback_P.KalibrerElev_Gain * rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopterExercise4_feedback_B.Add =
      helikopterExercise4_feedback_B.KalibrerElev +
      helikopterExercise4_feedback_P.Constant_Value;

    /* ToFile: '<Root>/To File1' */
    if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
      if (!(++helikopterExercise4_feedb_DWork.ToFile1_IWORK.Decimation % 1) &&
          (helikopterExercise4_feedb_DWork.ToFile1_IWORK.Count*2)+1 < 100000000 )
      {
        FILE *fp = (FILE *)
          helikopterExercise4_feedb_DWork.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4_feedb_DWork.ToFile1_IWORK.Decimation = 0;
          u[0] = helikopterExercise4_feedback_M->Timing.t[1];
          u[1] = helikopterExercise4_feedback_B.Add;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4_feedback_M,
                              "Error writing to MAT-file elevation.mat");
            return;
          }

          if (((++helikopterExercise4_feedb_DWork.ToFile1_IWORK.Count)*2)+1 >=
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

    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterExercise4_feedback_B.KalibrerPitch =
      helikopterExercise4_feedback_P.KalibrerPitch_Gain * rtb_HILReadEncoder_o2;

    /* ToFile: '<Root>/To File2' */
    if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
      if (!(++helikopterExercise4_feedb_DWork.ToFile2_IWORK.Decimation % 1) &&
          (helikopterExercise4_feedb_DWork.ToFile2_IWORK.Count*2)+1 < 100000000 )
      {
        FILE *fp = (FILE *)
          helikopterExercise4_feedb_DWork.ToFile2_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          helikopterExercise4_feedb_DWork.ToFile2_IWORK.Decimation = 0;
          u[0] = helikopterExercise4_feedback_M->Timing.t[1];
          u[1] = helikopterExercise4_feedback_B.KalibrerPitch;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(helikopterExercise4_feedback_M,
                              "Error writing to MAT-file pitch.mat");
            return;
          }

          if (((++helikopterExercise4_feedb_DWork.ToFile2_IWORK.Count)*2)+1 >=
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
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopterExercise4_feedback_X.Integrator_CSTATE >=
      helikopterExercise4_feedback_P.Integrator_UpperSat ) {
    helikopterExercise4_feedback_X.Integrator_CSTATE =
      helikopterExercise4_feedback_P.Integrator_UpperSat;
  } else if (helikopterExercise4_feedback_X.Integrator_CSTATE <=
             helikopterExercise4_feedback_P.Integrator_LowerSat ) {
    helikopterExercise4_feedback_X.Integrator_CSTATE =
      helikopterExercise4_feedback_P.Integrator_LowerSat;
  }

  rtb_Saturation_g = helikopterExercise4_feedback_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopterExercise4_feedback_B.KalibrerVandring =
      helikopterExercise4_feedback_P.KalibrerVandring_Gain *
      rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopterExercise4_feedback_P.VandringDeriv_D*
    helikopterExercise4_feedback_B.KalibrerVandring;
  rtb_VandringDeriv += helikopterExercise4_feedback_P.VandringDeriv_C*
    helikopterExercise4_feedback_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopterExercise4_feedback_P.TransferFcn4_D*
    helikopterExercise4_feedback_B.KalibrerPitch;
  rtb_Gain2 += helikopterExercise4_feedback_P.TransferFcn4_C*
    helikopterExercise4_feedback_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopterExercise4_feedback_P.TransferFcn5_D*
    helikopterExercise4_feedback_B.KalibrerElev;
  rtb_Gain1 += helikopterExercise4_feedback_P.TransferFcn5_C*
    helikopterExercise4_feedback_X.TransferFcn5_CSTATE;

  /* SignalConversion: '<Root>/TmpSignal ConversionAtdeg->radInport1' */
  rtb_Sum1_p[0] = helikopterExercise4_feedback_B.VandringLavpass;
  rtb_Sum1_p[1] = rtb_VandringDeriv;
  rtb_Sum1_p[2] = helikopterExercise4_feedback_B.KalibrerPitch;
  rtb_Sum1_p[3] = rtb_Gain2;
  rtb_Sum1_p[4] = helikopterExercise4_feedback_B.Add;
  rtb_Sum1_p[5] = rtb_Gain1;

  /* Gain: '<Root>/deg->rad' */
  for (tmp = 0; tmp < 6; tmp++) {
    rtb_degrad[tmp] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_degrad[tmp] += helikopterExercise4_feedback_P.degrad_Gain[6 * tmp_1 +
        tmp] * rtb_Sum1_p[tmp_1];
    }
  }

  /* Gain: '<S1>/K_ed' */
  rtb_K_ed = helikopterExercise4_feedback_P.K_ed_Gain * rtb_degrad[5];

  /* FromWorkspace: '<Root>/u_star' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise4_feedb_DWork.u_star_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise4_feedb_DWork.u_star_PWORK.TimePtr;
    int_T currTimeIndex = helikopterExercise4_feedb_DWork.u_star_IWORK.PrevIndex;
    real_T t = helikopterExercise4_feedback_M->Timing.t[0];

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

    helikopterExercise4_feedb_DWork.u_star_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Sum[0] = pDataValues[currTimeIndex];
          pDataValues += 141;
          rtb_Sum[1] = pDataValues[currTimeIndex];
          pDataValues += 141;
        } else {
          rtb_Sum[0] = pDataValues[currTimeIndex + 1];
          pDataValues += 141;
          rtb_Sum[1] = pDataValues[currTimeIndex + 1];
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
        rtb_Sum[0] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum[1] = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* FromWorkspace: '<Root>/x_star' */
  {
    real_T *pDataValues = (real_T *)
      helikopterExercise4_feedb_DWork.x_star_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterExercise4_feedb_DWork.x_star_PWORK.TimePtr;
    int_T currTimeIndex = helikopterExercise4_feedb_DWork.x_star_IWORK.PrevIndex;
    real_T t = helikopterExercise4_feedback_M->Timing.t[0];

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

    helikopterExercise4_feedb_DWork.x_star_IWORK.PrevIndex = currTimeIndex;

    /* post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T i1;
            real_T *y0 = rtb_Sum1_p;
            for (i1=0; i1 < 6; i1++) {
              y0[i1] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T i1;
            real_T *y0 = rtb_Sum1_p;
            for (i1=0; i1 < 6; i1++) {
              y0[i1] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
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
          real_T *y0 = rtb_Sum1_p;
          for (i1=0; i1 < 6; i1++) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            y0[i1] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* Sum: '<S3>/Sum1' */
  for (tmp = 0; tmp < 6; tmp++) {
    rtb_Sum1_p[tmp] = rtb_degrad[tmp] - rtb_Sum1_p[tmp];
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/Gain'
   */
  for (tmp = 0; tmp < 2; tmp++) {
    tmp_0[tmp] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      tmp_0[tmp] += helikopterExercise4_feedback_P.Gain_Gain[(tmp_1 << 1) + tmp]
        * rtb_Sum1_p[tmp_1];
    }

    rtb_Sum[tmp] -= tmp_0[tmp];
  }

  /* Sum: '<S1>/Sum' */
  rtb_Gain1 = rtb_Sum[1] - rtb_degrad[4];

  /* Gain: '<S1>/K_ei' */
  helikopterExercise4_feedback_B.K_ei = helikopterExercise4_feedback_P.K_ei_Gain
    * rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_g = (helikopterExercise4_feedback_P.K_ep_Gain * rtb_Gain1 +
                      rtb_Saturation_g) + rtb_K_ed;

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_g = rt_SATURATE(rtb_Saturation_g,
    helikopterExercise4_feedback_P.Saturation_LowerSat,
    helikopterExercise4_feedback_P.Saturation_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
  }

  /* Sum: '<S4>/Sum' incorporates:
   *  Gain: '<S4>/K_pd'
   *  Gain: '<S4>/K_pp'
   *  Saturate: '<S4>/Saturation'
   *  Sum: '<S4>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Sum[0],
    helikopterExercise4_feedback_P.Saturation_LowerSat_l,
    helikopterExercise4_feedback_P.Saturation_UpperSat_k) - rtb_degrad[2]) *
    helikopterExercise4_feedback_P.K_pp_Gain -
    helikopterExercise4_feedback_P.K_pd_Gain * rtb_degrad[3];

  /* Gain: '<S5>/Gain2' incorporates:
   *  Sum: '<S5>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_g - rtb_Gain1) *
    helikopterExercise4_feedback_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopterExercise4_feedback_B.SatB = rt_SATURATE(rtb_Gain2,
    helikopterExercise4_feedback_P.SatB_LowerSat,
    helikopterExercise4_feedback_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
  }

  /* Gain: '<S5>/Gain1' incorporates:
   *  Sum: '<S5>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_g) *
    helikopterExercise4_feedback_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopterExercise4_feedback_B.Sat = rt_SATURATE(rtb_Gain1,
    helikopterExercise4_feedback_P.Sat_LowerSat,
    helikopterExercise4_feedback_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopterExercise4_feedback/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopterExercise4_feedb_DWork.HILWriteAnalog_Buffer[0] =
        helikopterExercise4_feedback_B.SatB;
      helikopterExercise4_feedb_DWork.HILWriteAnalog_Buffer[1] =
        helikopterExercise4_feedback_B.Sat;
      result = hil_write_analog
        (helikopterExercise4_feedb_DWork.HILInitialize_Card,
         helikopterExercise4_feedback_P.HILWriteAnalog_Channels, 2,
         &helikopterExercise4_feedb_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopterExercise4_feedback_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
    rt_ertODEUpdateContinuousStates(&helikopterExercise4_feedback_M->solverInfo);
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
  if (!(++helikopterExercise4_feedback_M->Timing.clockTick0)) {
    ++helikopterExercise4_feedback_M->Timing.clockTickH0;
  }

  helikopterExercise4_feedback_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopterExercise4_feedback_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopterExercise4_feedback_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopterExercise4_feedback_M->Timing.clockTick1)) {
      ++helikopterExercise4_feedback_M->Timing.clockTickH1;
    }

    helikopterExercise4_feedback_M->Timing.t[1] =
      helikopterExercise4_feedback_M->Timing.clockTick1 *
      helikopterExercise4_feedback_M->Timing.stepSize1 +
      helikopterExercise4_feedback_M->Timing.clockTickH1 *
      helikopterExercise4_feedback_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopterExercise4_feedback_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->VandringLavpass_CSTATE =
      helikopterExercise4_feedback_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->VandringLavpass_CSTATE +=
      (helikopterExercise4_feedback_P.VandringLavpass_A)*
      helikopterExercise4_feedback_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopterExercise4_feedback_X.Integrator_CSTATE <=
            helikopterExercise4_feedback_P.Integrator_LowerSat );
    usat = ( helikopterExercise4_feedback_X.Integrator_CSTATE >=
            helikopterExercise4_feedback_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopterExercise4_feedback_B.K_ei > 0)) ||
        (usat && (helikopterExercise4_feedback_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopterExer *)
        helikopterExercise4_feedback_M->ModelData.derivs)->Integrator_CSTATE =
        helikopterExercise4_feedback_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopterExer *)
        helikopterExercise4_feedback_M->ModelData.derivs)->Integrator_CSTATE =
        0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->VandringDeriv_CSTATE =
      helikopterExercise4_feedback_B.KalibrerVandring;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->VandringDeriv_CSTATE +=
      (helikopterExercise4_feedback_P.VandringDeriv_A)*
      helikopterExercise4_feedback_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->TransferFcn4_CSTATE =
      helikopterExercise4_feedback_B.KalibrerPitch;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->TransferFcn4_CSTATE +=
      (helikopterExercise4_feedback_P.TransferFcn4_A)*
      helikopterExercise4_feedback_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->TransferFcn5_CSTATE =
      helikopterExercise4_feedback_B.KalibrerElev;
    ((StateDerivatives_helikopterExer *)
      helikopterExercise4_feedback_M->ModelData.derivs)->TransferFcn5_CSTATE +=
      (helikopterExercise4_feedback_P.TransferFcn5_A)*
      helikopterExercise4_feedback_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopterExercise4_feedback_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopterExercise4_feedback_P.Integrator_UpperSat = rtInf;
  helikopterExercise4_feedback_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopterExercise4_feedback_M, 0,
                sizeof(RT_MODEL_helikopterExercise4_fe));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopterExercise4_feedback_M->solverInfo,
                          &helikopterExercise4_feedback_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopterExercise4_feedback_M->solverInfo, &rtmGetTPtr
                (helikopterExercise4_feedback_M));
    rtsiSetStepSizePtr(&helikopterExercise4_feedback_M->solverInfo,
                       &helikopterExercise4_feedback_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopterExercise4_feedback_M->solverInfo,
                 &helikopterExercise4_feedback_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopterExercise4_feedback_M->solverInfo,
                         &helikopterExercise4_feedback_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopterExercise4_feedback_M->solverInfo,
      &helikopterExercise4_feedback_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopterExercise4_feedback_M->solverInfo,
                          (&rtmGetErrorStatus(helikopterExercise4_feedback_M)));
    rtsiSetRTModelPtr(&helikopterExercise4_feedback_M->solverInfo,
                      helikopterExercise4_feedback_M);
  }

  rtsiSetSimTimeStep(&helikopterExercise4_feedback_M->solverInfo,
                     MAJOR_TIME_STEP);
  helikopterExercise4_feedback_M->ModelData.intgData.f[0] =
    helikopterExercise4_feedback_M->ModelData.odeF[0];
  helikopterExercise4_feedback_M->ModelData.contStates = ((real_T *)
    &helikopterExercise4_feedback_X);
  rtsiSetSolverData(&helikopterExercise4_feedback_M->solverInfo, (void *)
                    &helikopterExercise4_feedback_M->ModelData.intgData);
  rtsiSetSolverName(&helikopterExercise4_feedback_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      helikopterExercise4_feedback_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopterExercise4_feedback_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopterExercise4_feedback_M->Timing.sampleTimes =
      (&helikopterExercise4_feedback_M->Timing.sampleTimesArray[0]);
    helikopterExercise4_feedback_M->Timing.offsetTimes =
      (&helikopterExercise4_feedback_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopterExercise4_feedback_M->Timing.sampleTimes[0] = (0.0);
    helikopterExercise4_feedback_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopterExercise4_feedback_M->Timing.offsetTimes[0] = (0.0);
    helikopterExercise4_feedback_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopterExercise4_feedback_M,
             &helikopterExercise4_feedback_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopterExercise4_feedback_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopterExercise4_feedback_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopterExercise4_feedback_M, -1);
  helikopterExercise4_feedback_M->Timing.stepSize0 = 0.001;
  helikopterExercise4_feedback_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopterExercise4_feedback_M->Sizes.checksums[0] = (21000766U);
  helikopterExercise4_feedback_M->Sizes.checksums[1] = (3298841055U);
  helikopterExercise4_feedback_M->Sizes.checksums[2] = (3763219128U);
  helikopterExercise4_feedback_M->Sizes.checksums[3] = (1228058123U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopterExercise4_feedback_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopterExercise4_feedback_M->extModeInfo,
      &helikopterExercise4_feedback_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopterExercise4_feedback_M->extModeInfo,
                        helikopterExercise4_feedback_M->Sizes.checksums);
    rteiSetTPtr(helikopterExercise4_feedback_M->extModeInfo, rtmGetTPtr
                (helikopterExercise4_feedback_M));
  }

  helikopterExercise4_feedback_M->solverInfoPtr =
    (&helikopterExercise4_feedback_M->solverInfo);
  helikopterExercise4_feedback_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopterExercise4_feedback_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopterExercise4_feedback_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopterExercise4_feedback_M->ModelData.blockIO = ((void *)
    &helikopterExercise4_feedback_B);

  {
    helikopterExercise4_feedback_B.VandringLavpass = 0.0;
    helikopterExercise4_feedback_B.KalibrerElev = 0.0;
    helikopterExercise4_feedback_B.Add = 0.0;
    helikopterExercise4_feedback_B.KalibrerPitch = 0.0;
    helikopterExercise4_feedback_B.KalibrerVandring = 0.0;
    helikopterExercise4_feedback_B.K_ei = 0.0;
    helikopterExercise4_feedback_B.SatB = 0.0;
    helikopterExercise4_feedback_B.Sat = 0.0;
  }

  /* parameters */
  helikopterExercise4_feedback_M->ModelData.defaultParam = ((real_T *)
    &helikopterExercise4_feedback_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopterExercise4_feedback_X;
    helikopterExercise4_feedback_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopterExercise4_feedback_X, 0,
                  sizeof(ContinuousStates_helikopterExer));
  }

  /* states (dwork) */
  helikopterExercise4_feedback_M->Work.dwork = ((void *)
    &helikopterExercise4_feedb_DWork);
  (void) memset((void *)&helikopterExercise4_feedb_DWork, 0,
                sizeof(D_Work_helikopterExercise4_feed));

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise4_feedb_DWork.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterExercise4_feedb_DWork.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopterExercise4_feedb_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopterExercise4_feedb_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopterExercise4_feedback_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopterExercise4_feedback_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise4_feedback/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopterExercise4_feedb_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopterExercise4_feedb_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopterExercise4_feedb_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopterExercise4_feedb_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopterExercise4_feedback_P.HILInitialize_AOTerminate &&
         !is_switching) || (helikopterExercise4_feedback_P.HILInitialize_AOExit &&
         is_switching)) {
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise4_feedback_P.HILInitialize_AOFinal;
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise4_feedback_P.HILInitialize_AOFinal;
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise4_feedback_P.HILInitialize_AOFinal;
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise4_feedback_P.HILInitialize_AOFinal;
      result = hil_write_analog
        (helikopterExercise4_feedb_DWork.HILInitialize_Card,
         &helikopterExercise4_feedback_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
      }
    }

    hil_close(helikopterExercise4_feedb_DWork.HILInitialize_Card);
    helikopterExercise4_feedb_DWork.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helikopterExercise4_feedb_DWork.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4_feedb_DWork.ToFile_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error writing header for ans to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error closing MAT-file travel.mat");
        return;
      }

      helikopterExercise4_feedb_DWork.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) helikopterExercise4_feedb_DWork.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "elevation.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error closing MAT-file elevation.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error reopening MAT-file elevation.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4_feedb_DWork.ToFile1_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error writing header for ans to MAT-file elevation.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error closing MAT-file elevation.mat");
        return;
      }

      helikopterExercise4_feedb_DWork.ToFile1_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File2' */
  {
    FILE *fp = (FILE *) helikopterExercise4_feedb_DWork.ToFile2_PWORK.FilePtr;
    if (fp != (NULL)) {
      const char *fileName = "pitch.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error reopening MAT-file pitch.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2,
           helikopterExercise4_feedb_DWork.ToFile2_IWORK.Count, "ans")) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error writing header for ans to MAT-file pitch.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helikopterExercise4_feedback_M,
                          "Error closing MAT-file pitch.mat");
        return;
      }

      helikopterExercise4_feedb_DWork.ToFile2_PWORK.FilePtr = (NULL);
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
  helikopterExercise4_feedback_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopterExercise4_feedback_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopterExercise4_feedback_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopterExercise4_feedback_M->Sizes.numY = (0);/* Number of model outputs */
  helikopterExercise4_feedback_M->Sizes.numU = (0);/* Number of model inputs */
  helikopterExercise4_feedback_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopterExercise4_feedback_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopterExercise4_feedback_M->Sizes.numBlocks = (49);/* Number of blocks */
  helikopterExercise4_feedback_M->Sizes.numBlockIO = (8);/* Number of block outputs */
  helikopterExercise4_feedback_M->Sizes.numBlockPrms = (160);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterExercise4_feedback_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopterExercise4_feedback_X.Integrator_CSTATE =
    helikopterExercise4_feedback_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopterExercise4_feedback_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopterExercise4_feedback_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopterExercise4_feedback_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterExercise4_feedback/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("sensoray_model_626", "0",
                      &helikopterExercise4_feedb_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopterExercise4_feedback_P.HILInitialize_CKPStart && !is_switching)
        || (helikopterExercise4_feedback_P.HILInitialize_CKPEnter &&
            is_switching)) {
      {
        int_T i1;
        int32_T *dw_ClockModes =
          &helikopterExercise4_feedb_DWork.HILInitialize_ClockModes[0];
        for (i1=0; i1 < 6; i1++) {
          dw_ClockModes[i1] =
            helikopterExercise4_feedback_P.HILInitialize_CKModes;
        }
      }

      result = hil_set_clock_mode
        (helikopterExercise4_feedb_DWork.HILInitialize_Card, (t_clock *)
         &helikopterExercise4_feedback_P.HILInitialize_CKChannels[0], 6U,
         (t_clock_mode *)
         &helikopterExercise4_feedb_DWork.HILInitialize_ClockModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear
      (helikopterExercise4_feedb_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
      return;
    }

    if ((helikopterExercise4_feedback_P.HILInitialize_AIPStart && !is_switching)
        || (helikopterExercise4_feedback_P.HILInitialize_AIPEnter &&
            is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helikopterExercise4_feedb_DWork.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMinimums[i1] = helikopterExercise4_feedback_P.HILInitialize_AILow;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helikopterExercise4_feedb_DWork.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMaximums[i1] =
            helikopterExercise4_feedback_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helikopterExercise4_feedb_DWork.HILInitialize_Card,
         &helikopterExercise4_feedback_P.HILInitialize_AIChannels[0], 16U,
         &helikopterExercise4_feedb_DWork.HILInitialize_AIMinimums[0],
         &helikopterExercise4_feedb_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise4_feedback_P.HILInitialize_AOPStart && !is_switching)
        || (helikopterExercise4_feedback_P.HILInitialize_AOPEnter &&
            is_switching)) {
      helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[0] =
        helikopterExercise4_feedback_P.HILInitialize_AOLow;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[1] =
        helikopterExercise4_feedback_P.HILInitialize_AOLow;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[2] =
        helikopterExercise4_feedback_P.HILInitialize_AOLow;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[3] =
        helikopterExercise4_feedback_P.HILInitialize_AOLow;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[0] =
        helikopterExercise4_feedback_P.HILInitialize_AOHigh;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[1] =
        helikopterExercise4_feedback_P.HILInitialize_AOHigh;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[2] =
        helikopterExercise4_feedback_P.HILInitialize_AOHigh;
      helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[3] =
        helikopterExercise4_feedback_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopterExercise4_feedb_DWork.HILInitialize_Card,
         &helikopterExercise4_feedback_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4_feedb_DWork.HILInitialize_AOMinimums[0],
         &helikopterExercise4_feedb_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterExercise4_feedback_P.HILInitialize_AOStart && !is_switching) ||
        (helikopterExercise4_feedback_P.HILInitialize_AOEnter && is_switching))
    {
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[0] =
        helikopterExercise4_feedback_P.HILInitialize_AOInitial;
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[1] =
        helikopterExercise4_feedback_P.HILInitialize_AOInitial;
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[2] =
        helikopterExercise4_feedback_P.HILInitialize_AOInitial;
      helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[3] =
        helikopterExercise4_feedback_P.HILInitialize_AOInitial;
      result = hil_write_analog
        (helikopterExercise4_feedb_DWork.HILInitialize_Card,
         &helikopterExercise4_feedback_P.HILInitialize_AOChannels[0], 4U,
         &helikopterExercise4_feedb_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterExercise4_feedback_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    const char *fileName = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4_feedback_M,
                        "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4_feedback_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    helikopterExercise4_feedb_DWork.ToFile_IWORK.Count = 0;
    helikopterExercise4_feedb_DWork.ToFile_IWORK.Decimation = -1;
    helikopterExercise4_feedb_DWork.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File1' */
  {
    const char *fileName = "elevation.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4_feedback_M,
                        "Error creating .mat file elevation.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4_feedback_M,
                        "Error writing mat file header to file elevation.mat");
      return;
    }

    helikopterExercise4_feedb_DWork.ToFile1_IWORK.Count = 0;
    helikopterExercise4_feedb_DWork.ToFile1_IWORK.Decimation = -1;
    helikopterExercise4_feedb_DWork.ToFile1_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File2' */
  {
    const char *fileName = "pitch.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helikopterExercise4_feedback_M,
                        "Error creating .mat file pitch.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"ans")) {
      rtmSetErrorStatus(helikopterExercise4_feedback_M,
                        "Error writing mat file header to file pitch.mat");
      return;
    }

    helikopterExercise4_feedb_DWork.ToFile2_IWORK.Count = 0;
    helikopterExercise4_feedb_DWork.ToFile2_IWORK.Decimation = -1;
    helikopterExercise4_feedb_DWork.ToFile2_PWORK.FilePtr = fp;
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
      1.6455121993179131E-001, 8.2275609965895669E-002, 3.7536039988674536E-001,
      4.3316674352228418E-001, 4.3750730348445968E-002, -6.3032740765410028E-017,
      2.9320995138632851E-017, -4.9242001927076724E-017, 2.7492193328767314E-017,
      -2.4243789811287732E-017, 1.2547694178817646E-017, 4.3140054752329482E-017,
      1.8333696404830535E-017, 2.0670836564523458E-017, 7.5131971967406885E-017,
      1.5248343219230958E-017, 2.5892519912054760E-017, 5.2616878258057187E-017,
      -2.8156348487569725E-017, 5.1479873722659928E-017,
      -3.6313514575414520E-017, 6.4764275503546270E-017, 2.5274634354805459E-017,
      4.2060498540172030E-017, 6.6205561288749272E-017, -3.7789495398254590E-017,
      2.6618321866493838E-017, -5.8611677465715666E-018,
      -7.1531039932645045E-017, 7.3666410712737186E-017,
      -1.8553782067976785E-017, 4.5062555520040899E-017, 5.4752898633100480E-018,
      3.0127820000582697E-017, -8.1380728448693401E-017,
      -1.4848648183569015E-017, -2.2142179465426977E-017,
      -8.5824031720201439E-017, 6.7385368066366118E-017,
      -3.5472912767166068E-017, -8.3999406179079832E-018,
      1.2301874037014872E-018, -3.6853485497152337E-017, 3.1504123257711109E-017,
      -2.9548340367001838E-018, 2.0129902996921763E-017,
      -2.1489072704774011E-017, 2.4254564017166306E-017, 3.7474079912795470E-017,
      5.2718269552819312E-017, 5.3591167772723006E-017, 2.7473500875308096E-017,
      4.9917538984345293E-017, 2.8069070287503149E-017, 2.8828644367399364E-017,
      -8.4084070830440811E-017, 8.6446915939149821E-017,
      -4.1886171853938351E-017, 6.9353600220097777E-017,
      -5.0253479098508122E-017, 4.1780685444203925E-017,
      -7.0200151674595011E-017, -1.3444441563470694E-017,
      -5.5110921477500673E-017, 2.2261657353571109E-017,
      -4.2029596037474950E-017, 3.3032927049418158E-017,
      -2.3663449853943877E-018, -2.3396878382724158E-017,
      -8.8335481503882100E-018, -4.2620158265303767E-017,
      -4.7798831101750747E-017, 5.1238680488941411E-017, 4.3703547467713110E-019,
      5.0538312395771987E-017, 3.7868725617080016E-018, 4.8736473726019412E-018,
      -2.0311781785150003E-017, -4.0252830811214003E-017,
      -2.9932962735486687E-017, -3.2929544620522013E-017,
      -4.0314541948298291E-017, -6.9566735276824713E-017,
      -3.6959913026942634E-017, -3.4908059156007092E-017,
      -2.0640556154797819E-017, 7.6293080915830684E-017,
      -1.8261918955217197E-017, 2.7224500461921474E-017,
      -1.1446868264773282E-017, 3.0160218979976284E-017,
      -2.8484368361797763E-017, 5.0915674721959847E-017, 4.5507347003341675E-017,
      4.8185499349935732E-018, 2.0757101868887950E-017, -1.3469698983735974E-017,
      -7.2982996310129084E-005, -1.1763922509238791E-005,
      -1.9288246651854958E-004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 2.1688164126099841E-006, 4.2100908691366903E-006,
      4.2440041831924772E-006, 8.7792050716119898E-006, 1.0084162191669228E-005,
      1.4637647505480010E-005, 1.8206034407835506E-005, 2.6592961709895451E-005,
      3.6607033005695992E-005, 5.0702482063060162E-005, 7.3613851237233074E-005,
      9.7731745477375485E-005, 1.3764168551811851E-004, 1.9173158630859610E-004,
      2.6742002791671938E-004, 3.6960093737397573E-004, 5.1463897680764940E-004,
      7.1475030475876130E-004, 9.9552206361313775E-004, 1.3836895161814017E-003,
      1.9246948502176383E-003, 2.6733444718873801E-003, 3.7136894076575042E-003,
      5.1535692301117163E-003, 7.1442618510898590E-003, 9.8806949037807540E-003,
      1.3632115365728196E-002, 1.8741524551384534E-002, 2.5640607704581118E-002,
      3.4855265347152199E-002, 4.6957223554381457E-002, 6.2484048706435411E-002,
      8.1620067507432392E-002, 1.0361117376851657E-001, 1.2940269742213431E-001,
      1.6164280929722993E-001, 1.9039389648954430E-001, 2.0028084672216104E-001,
      1.5951186212424165E-001, -1.6817045889495980E-004, 6.8288064581008192E-005,
      -1.9153128611288945E-004, 1.4969755201495617E-005, 1.0301307251793117E-004,
      7.4834194349023063E-005, 2.1965973151294178E-005, -8.4370319420995195E-006,
      -1.6869014401749423E-005, -1.3574233390748299E-005,
      -8.4167876672468109E-006, -2.9559012936429271E-006,
      1.6083082402109321E-006, 3.5265656923943928E-006, 2.8196924945167159E-006,
      3.5443203721843075E-006, 4.2655376887169274E-006, 3.5590380075127053E-006,
      2.5154977880514868E-006, 2.7359817593852697E-006, 1.2063758063773728E-006,
      1.4056292202521963E-006, 2.2255777697058449E-006, 5.1071720091987753E-007,
      6.4691486020812893E-007, 1.8002985481944800E-006, 3.5876915501074068E-007,
      1.5872483368133255E-006, -2.3096780013119976E-007,
      -2.1216450551048701E-007, 9.2265119339001373E-008,
      -8.5641840965893249E-007, 1.2288176761401515E-006, 7.6925071259704750E-007,
      -1.0322416066339622E-006, 1.3417369148199057E-006, 1.1424713947549451E-006,
      -1.1947067506621205E-007, -1.7506030270237668E-007,
      1.3221612566856921E-006, -2.6148130517918231E-007, 8.7255312256192961E-007,
      7.6440902378459122E-007, 7.2573718102500926E-007, 7.1190432020145378E-007,
      7.0676581625651423E-007, 7.0446830784805734E-007, 7.0286285501261138E-007,
      7.0112257152089517E-007, 6.9888886403394108E-007, 6.9599227117109262E-007,
      6.9244367862840143E-007, 6.8870012210321063E-007, 6.8652371103632696E-007,
      6.9142613058514114E-007, 7.1948304371324691E-007, 8.1622774853322162E-007,
      1.1089298747294580E-006, -1.2343101349720172E-006,
      -1.0868553586475565E-006, -3.4128195669609266E-019, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helikopterExercise4_feedb_DWork.u_star_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise4_feedb_DWork.u_star_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise4_feedb_DWork.u_star_IWORK.PrevIndex = 0;
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
      3.1415926535897931E+000, 3.1415926535897931E+000, 3.1296829880716439E+000,
      3.1058636570353455E+000, 3.0808636570353456E+000, 3.0558636570353457E+000,
      3.0308636570353458E+000, 3.0058636570353459E+000, 2.9808636570353459E+000,
      2.9558636570353460E+000, 2.9308636570353461E+000, 2.9058636570353462E+000,
      2.8808636570353463E+000, 2.8558636570353464E+000, 2.8308636570353465E+000,
      2.8058636570353466E+000, 2.7808636570353467E+000, 2.7558636570353467E+000,
      2.7308636570353468E+000, 2.7058636570353469E+000, 2.6808636570353470E+000,
      2.6558636570353471E+000, 2.6308636570353472E+000, 2.6058636570353473E+000,
      2.5808636570353474E+000, 2.5558636570353475E+000, 2.5308636570353475E+000,
      2.5058636570353476E+000, 2.4808636570353477E+000, 2.4558636570353478E+000,
      2.4308636570353479E+000, 2.4058636570353480E+000, 2.3808636570353481E+000,
      2.3558636570353482E+000, 2.3308636570353483E+000, 2.3058636570353483E+000,
      2.2808636570353484E+000, 2.2558636570353485E+000, 2.2308636570353486E+000,
      2.2058636570353487E+000, 2.1808636570353488E+000, 2.1558636570353489E+000,
      2.1308636570353490E+000, 2.1058636570353491E+000, 2.0808636570353491E+000,
      2.0558636570353492E+000, 2.0308636570353493E+000, 2.0058636570353494E+000,
      1.9808636570353495E+000, 1.9558636570353496E+000, 1.9308636570353497E+000,
      1.9058636570353498E+000, 1.8808636570353499E+000, 1.8558636570353499E+000,
      1.8308636570353500E+000, 1.8058636570353501E+000, 1.7808636570353502E+000,
      1.7558636570353503E+000, 1.7308636570353504E+000, 1.7058636570353505E+000,
      1.6808636570353506E+000, 1.6558636570353507E+000, 1.6308636570353507E+000,
      1.6058636570353508E+000, 1.5808636570353509E+000, 1.5558636570353510E+000,
      1.5308636570353511E+000, 1.5058636570353512E+000, 1.4808636570353513E+000,
      1.4558636570353514E+000, 1.4308636570353515E+000, 1.4058636570353515E+000,
      1.3808636570353516E+000, 1.3558636570353517E+000, 1.3308636570353518E+000,
      1.3058636570353519E+000, 1.2808636570353520E+000, 1.2558636570353521E+000,
      1.2308636570353522E+000, 1.2058636570353523E+000, 1.1808636570353523E+000,
      1.1558636570353524E+000, 1.1308636570353525E+000, 1.1058636570353526E+000,
      1.0808636570353527E+000, 1.0558636570353528E+000, 1.0308636570353529E+000,
      1.0058636570353530E+000, 9.8086365703535294E-001, 9.5586365703535292E-001,
      9.3086365703535290E-001, 9.0586365703535288E-001, 8.8086365703535285E-001,
      8.5586365703535283E-001, 8.3086365703535281E-001, 8.0586365703535279E-001,
      7.8086365703535277E-001, 7.5586365703535274E-001, 7.3086365703535272E-001,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -4.7638662072596753E-002, -9.5277324145193507E-002, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
      -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -9.9978870941827683E-002, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.1907152458564879E-002,
      -5.5739886787617387E-017, 2.0772623571294503E-017,
      -4.0689350866713522E-017, 1.7242388117443004E-017,
      -3.8937850403208345E-017, 2.1766949325576189E-017, 3.8170763955892774E-017,
      6.3336639780572349E-019, 3.5214251666166543E-017, 7.4245914229417526E-017,
      8.4532194678700323E-018, 4.4168985306391686E-017, 4.7933053763281726E-017,
      -1.2104337239219392E-017, 3.2818196541676402E-017,
      -4.0720257373032180E-017, 5.5430062757766787E-017, 2.1377072410499187E-017,
      6.8104850667236576E-017, 5.7142927702979019E-017, -2.0248731817851526E-017,
      4.7284221430145724E-017, -4.2296623130871562E-017,
      -6.3264690593865556E-017, 6.3156252237994379E-017,
      -2.6357581706921172E-017, 5.1842454745298165E-017, 3.7573974915922674E-017,
      2.3134707599876680E-017, -7.1783049147531808E-017,
      -2.0722233732022353E-018, -6.7467264789140119E-017,
      -7.5424078883457990E-017, 7.2777866519854880E-017,
      -3.4415553468768622E-017, 2.0223315175436576E-018,
      -1.9561107660826181E-017, -4.1932854932213813E-017,
      3.5410514487480618E-017, 1.2901284350209852E-017, 1.1725838640425677E-017,
      -3.8140298260745783E-017, 1.3534339221869388E-017, 4.0680895756186534E-017,
      6.1219163907721745E-017, 4.0746004302001894E-017, 2.3932903904603018E-017,
      7.3212513983487299E-017, 5.9838819517498656E-017, -7.5984898020968318E-018,
      -7.4970558169365758E-017, 7.4625907171262027E-017,
      -2.9171345212892228E-017, 6.0341500550260852E-017,
      -2.4873252959043582E-017, 4.6410449347326967E-017,
      -6.9990786856543821E-017, -2.6521987284967919E-017,
      -4.8922581842392304E-017, 3.8902111105850257E-018,
      -4.2207361507121379E-017, 3.9299008785519946E-017,
      -8.8853010008558231E-019, -7.3331064908320332E-019,
      -5.5132985449156358E-018, -7.5238335533521120E-017,
      -5.0940786589895650E-017, 4.0700843451647878E-017, 4.8776695742170795E-018,
      5.6320448639423064E-017, 1.6469047926331949E-017, 1.7996103418541687E-017,
      -2.3121859229914750E-017, -4.4318277357265339E-017,
      -2.4778456617753713E-017, -2.1612195291820030E-017,
      -3.9596863987611255E-017, -6.6616114688935525E-017,
      -5.5390126625468086E-017, -7.0391937876024668E-017,
      -4.2464516866939098E-018, 7.5589815914405279E-017,
      -1.6494038870090617E-017, 3.1785664053813143E-017, 2.0204526075318262E-018,
      2.3895246539974299E-018, -4.0647763978734707E-017, 6.9236688580344304E-017,
      5.2998150154263573E-017, 2.8203102486603119E-018, 4.0736405309849764E-018,
      9.8053351367131243E-018, 6.0025936739475573E-017, -3.9145069127035878E-017,
      -2.3223047221010329E-004, -1.5354780831330797E-004, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0943951023931953E+000,
      1.8133581149616251E-016, -1.8867664925589358E+000,
      -2.0762860983425974E-001, 3.0605004143564747E-016,
      -2.4584789775203211E-016, 2.3172695593662612E-016,
      -2.2472095408260541E-016, 2.4281919891513814E-016, 6.5615258521266314E-017,
      -1.5014959023234822E-016, 1.3832354107344327E-016, 1.5612665025300391E-016,
      -2.6317077904618996E-016, 1.4286306335408660E-016, 1.5056273827560158E-017,
      -2.4014956401000446E-016, 1.7969013512358316E-016,
      -2.9415381565883430E-016, 3.8460128052319587E-016,
      -1.3621196138907040E-016, 1.8691111302694956E-016,
      -4.3847691857030220E-017, -3.0956663808332219E-016,
      2.7013181299198896E-016, -3.5832337824406914E-016,
      -8.3872269851975964E-017, 5.0568377132743979E-016,
      -3.5805533577966219E-016, 3.1280014580887736E-016,
      -5.7073919317501955E-017, -5.7757069264183986E-017,
      -3.7967102698963389E-016, 2.7884330309731829E-016,
      -2.6158016566375153E-016, -3.1827256377271452E-017,
      5.9280778161325148E-016, -4.2877367995449401E-016, 1.4575153994524913E-016,
      -8.6333756713479353E-017, -8.9486989085550516E-017,
      3.0937347767877772E-016, -9.0036920549083063E-017,
      -4.7017828391367042E-018, -1.9946454760468585E-016,
      2.0669854993046070E-016, 1.0858622613726857E-016, 8.2153072606140892E-017,
      -8.1892638422879441E-017, -6.7252401589595506E-017,
      1.9711844031553714E-016, -5.3494777863954566E-017,
      -2.6974923727838199E-016, -2.6948827346907572E-016,
      5.9838586136251109E-016, -4.1518900953661697E-016, 3.5805138305261232E-016,
      -3.4085901403721771E-016, 2.8513480922548220E-016,
      -4.6560494481548318E-016, 1.7387519828630362E-016,
      -8.9602378229697553E-017, 2.1125117181190934E-016,
      -1.8439029047082564E-016, 3.2602548117056533E-016,
      -1.6075015554242213E-016, 6.2087780400951661E-019,
      -1.9119951583329730E-017, -2.7890014795442192E-016,
      9.7190195774501871E-017, 3.6656652016617414E-016, -1.4329269550972320E-016,
      2.0577111626082395E-016, -1.5940560285236446E-016, 6.1082219688389508E-018,
      -1.6447185059382575E-016, -8.4785672509402367E-017,
      7.8159282958046516E-017, 1.2665045303734732E-017, -7.1938674783164911E-017,
      -1.0807700280529708E-016, 4.4903952253869737E-017,
      -6.0007245002226341E-017, 2.6458194475732306E-016, 3.1934507040439673E-016,
      -3.6833541913798356E-016, 1.9311881169561504E-016,
      -1.1906084578512526E-016, 1.4762881858624150E-018,
      -1.7214915453092854E-016, 4.3953781023631599E-016,
      -6.4954153704322887E-017, -2.0071135962241305E-016,
      5.0133211292986579E-018, 2.2926778422912595E-017, 2.0088240641104983E-016,
      -3.9668402346604580E-016, -9.2892188884025647E-004,
      3.1473065558718127E-004, 3.4345302191836975E-004, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.8189390424477147E-007,
      9.0163012564646284E-007, 1.5760550718110372E-006, 2.7734664732707130E-006,
      4.1873860730894379E-006, 6.0932177913065728E-006, 8.4056377552765693E-006,
      1.1665051281866372E-005, 1.6169136957996390E-005, 2.2401876133643521E-005,
      3.1471867662354919E-005, 4.3596468254030887E-005, 6.0515410715914030E-005,
      8.4122322118424848E-005, 1.1708853690610523E-004, 1.6267544518280968E-004,
      2.2607604968514410E-004, 3.1414422360102494E-004, 4.3681189666925859E-004,
      6.0738702432772520E-004, 8.4466162587234128E-004, 1.1742324623471971E-003,
      1.6319291213071896E-003, 2.2669034816359523E-003, 3.1467404895038856E-003,
      4.3627148967588869E-003, 6.0384075872087209E-003, 8.3384295724903256E-003,
      1.1477994261789377E-002, 1.5732293008596805E-002, 2.1438286381073330E-002,
      2.8982929394823866E-002, 3.8746200933712407E-002, 5.0957974581594345E-002,
      6.5882986781963346E-002, 8.4109263915647131E-002, 1.0498194882756924E-001,
      1.2545159110278151E-001, 1.3780560579681475E-001, 1.2465650865428550E-001,
      1.0337096846383818E-001, 8.1667376108027212E-002, 6.2649674884482240E-002,
      4.7155388841440486E-002, 3.5035789306209736E-002, 2.5791443153531046E-002,
      1.8858147955179254E-002, 1.3719887452015099E-002, 9.9450439666232751E-003,
      7.1895015199348721E-003, 5.1875538530991661E-003, 3.7382388007271752E-003,
      2.6915596393879818E-003, 1.9367550066512185E-003, 1.3931812584053870E-003,
      1.0021532570227377E-003, 7.2093196584953215E-004, 5.1865038487208651E-004,
      3.7325910260734989E-004, 2.6859750529427591E-004, 1.9333812393557124E-004,
      1.3935341922896610E-004, 1.0040128990272143E-004, 7.2351162581546559E-005,
      5.2318658315055337E-005, 3.7807408806911024E-005, 2.7480086151170475E-005,
      1.9878994820892417E-005, 1.4324056781503206E-005, 1.0323081362899903E-005,
      7.3206123339612287E-006, 5.3661004242613309E-006, 4.0117234435086167E-006,
      2.8316383201607561E-006, 2.1809915166371874E-006, 1.7940417804068708E-006,
      1.3954824857632624E-006, 1.0370051128020461E-006, 9.3524589229461467E-007,
      7.4029326815937673E-007, 6.8198721539518067E-007, 6.7012067430207787E-007,
      6.7275408699669428E-007, 6.7886248038239030E-007, 6.8485475642820479E-007,
      6.8972459294792240E-007, 6.9331812591433479E-007, 6.9572413448794545E-007,
      6.9706721530202679E-007, 6.9744761507456401E-007, 6.9694453608260968E-007,
      6.9567814025458663E-007, 6.9399786504077408E-007, 6.9301290418640250E-007,
      6.9607195249022529E-007, 7.1287394837847457E-007, 7.7085836120632022E-007,
      5.3268957931557101E-007, 2.2995535738272132E-007, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1275756169791757E-006,
      2.4789448856067189E-006, 2.6976997846583721E-006, 4.7896456058386117E-006,
      5.6556783992748334E-006, 7.6233268728684565E-006, 9.2496798558800232E-006,
      1.3037654106359145E-005, 1.8016342704520075E-005, 2.4930956702588790E-005,
      3.6279966114845672E-005, 4.8498402366703760E-005, 6.7675769847532705E-005,
      9.4427645610043219E-005, 1.3186485915072158E-004, 1.8234763310681743E-004,
      2.5360241800933745E-004, 3.5227269566352325E-004, 4.9067069227293452E-004,
      6.8230051063386673E-004, 9.4909840617846478E-004, 1.3182833458994239E-003,
      1.8307866358399701E-003, 2.5398974413150502E-003, 3.5193480314717335E-003,
      4.8638976290200043E-003, 6.7027707617993357E-003, 9.2000879411264155E-003,
      1.2558258757196208E-002, 1.7017194987229714E-002, 2.2823973489906099E-002,
      3.0178572055002160E-002, 3.9053086155554154E-002, 4.8847094591527758E-002,
      5.9700048801476002E-002, 7.2905108534735197E-002, 8.3490739647688383E-002,
      8.1878569100849036E-002, 4.9416058776132983E-002, -5.2596388570117010E-002,
      -8.5142160761789307E-002, -8.6814369423243845E-002,
      -7.6070804894179958E-002, -6.1977144172166994E-002,
      -4.8478398140923001E-002, -3.6977384610714772E-002,
      -2.7733180793407176E-002, -2.0553042012656626E-002,
      -1.5099373941567290E-002, -1.1022169786753615E-002,
      -8.0077906673428258E-003, -5.7972602094879663E-003,
      -4.1867166453567733E-003, -3.0192185309470538E-003,
      -2.1742949929833260E-003, -1.5641120055305973E-003,
      -1.1248851646928220E-003, -8.0912632390978242E-004,
      -5.8156512905894661E-004, -4.1864638925229582E-004,
      -3.0103752543481869E-004, -2.1593881882642064E-004,
      -1.5580851730497857E-004, -1.1220050928469944E-004,
      -8.0130017065964900E-005, -5.8044998032577326E-005,
      -4.1309290622962308E-005, -3.0404365321112221E-005,
      -2.2219752157556845E-005, -1.6003901674413223E-005,
      -1.2009876115754629E-005, -7.8180476387995438E-006,
      -5.4175079230108414E-006, -4.7203404933913900E-006,
      -2.6025872140943208E-006, -1.5477989449212815E-006,
      -1.5942371785744538E-006, -1.4339094918449056E-006,
      -4.0703688202973111E-007, -7.7981049654098797E-007,
      -2.3322421105679774E-007, -4.7466164372393743E-008,
      1.0533650778516486E-008, 2.4433573542785653E-008, 2.3969104183137895E-008,
      1.9479346078787798E-008, 1.4374131865608637E-008, 9.6240342944290478E-009,
      5.3723232563097450E-009, 1.5215990902819661E-009, -2.0123159678326399E-009,
      -5.0655833121843125E-009, -6.7211008552634769E-009,
      -3.9398434174553247E-009, 1.2236193215290963E-008, 6.7207983553031907E-008,
      2.3193765131145608E-007, -9.5267512756292569E-007,
      -1.2109368877313170E-006, -5.8849767811499436E-007, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 } ;

    helikopterExercise4_feedb_DWork.x_star_PWORK.TimePtr = (void *) pTimeValues;
    helikopterExercise4_feedb_DWork.x_star_PWORK.DataPtr = (void *) pDataValues;
    helikopterExercise4_feedb_DWork.x_star_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopterExercise4_feedback_terminate();
}

RT_MODEL_helikopterExercise4_fe *helikopterExercise4_feedback(void)
{
  helikopterExercise4_feedback_initialize(1);
  return helikopterExercise4_feedback_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
