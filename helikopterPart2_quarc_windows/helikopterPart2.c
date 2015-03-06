/*
 * helikopterPart2.c
 *
 * Real-Time Workshop code generation for Simulink model "helikopterPart2.mdl".
 *
 * Model version              : 1.60
 * Real-Time Workshop version : 7.5  (R2010a)  25-Jan-2010
 * C source code generated on : Fri Mar 06 15:34:09 2015
 *
 * Target selection: quarc_windows.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helikopterPart2.h"
#include "helikopterPart2_private.h"
#include <stdio.h>
#include "helikopterPart2_dt.h"

/* Block signals (auto storage) */
BlockIO_helikopterPart2 helikopterPart2_B;

/* Continuous states */
ContinuousStates_helikopterPart helikopterPart2_X;

/* Block states (auto storage) */
D_Work_helikopterPart2 helikopterPart2_DWork;

/* Real-time model */
RT_MODEL_helikopterPart2 helikopterPart2_M_;
RT_MODEL_helikopterPart2 *helikopterPart2_M = &helikopterPart2_M_;

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
  helikopterPart2_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopterPart2_output(int_T tid)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T rtb_HILReadEncoder_o3;
  real_T rtb_VandringDeriv;
  real_T rtb_Gain2;
  real_T rtb_Saturation_h;
  real_T rtb_Gain1;
  real_T rtb_Gain[6];
  real_T tmp[6];
  int32_T tmp_0;
  int32_T tmp_1;
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    /* set solver stop time */
    if (!(helikopterPart2_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopterPart2_M->solverInfo,
                            ((helikopterPart2_M->Timing.clockTickH0 + 1) *
        helikopterPart2_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopterPart2_M->solverInfo,
                            ((helikopterPart2_M->Timing.clockTick0 + 1) *
        helikopterPart2_M->Timing.stepSize0 +
        helikopterPart2_M->Timing.clockTickH0 *
        helikopterPart2_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopterPart2_M)) {
    helikopterPart2_M->Timing.t[0] = rtsiGetT(&helikopterPart2_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    /* S-Function (hil_read_encoder_block): '<S2>/HIL Read Encoder' */

    /* S-Function Block: helikopterPart2/Heli 3D/HIL Read Encoder (hil_read_encoder_block) */
    {
      t_error result = hil_read_encoder(helikopterPart2_DWork.HILInitialize_Card,
        helikopterPart2_P.HILReadEncoder_Channels, 3,
        &helikopterPart2_DWork.HILReadEncoder_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
      } else {
        rtb_HILReadEncoder_o1 = helikopterPart2_DWork.HILReadEncoder_Buffer[0];
        rtb_HILReadEncoder_o2 = helikopterPart2_DWork.HILReadEncoder_Buffer[1];
        rtb_HILReadEncoder_o3 = helikopterPart2_DWork.HILReadEncoder_Buffer[2];
      }
    }

    /* Gain: '<S2>/Kalibrer-Elev' */
    helikopterPart2_B.KalibrerElev = helikopterPart2_P.KalibrerElev_Gain *
      rtb_HILReadEncoder_o3;

    /* Sum: '<Root>/Add' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopterPart2_B.Add = helikopterPart2_B.KalibrerElev +
      helikopterPart2_P.Constant_Value;
  }

  /* TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterPart2_B.VandringLavpass = helikopterPart2_P.VandringLavpass_C*
    helikopterPart2_X.VandringLavpass_CSTATE;
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    /* Gain: '<S2>/Kalibrer-Pitch' */
    helikopterPart2_B.KalibrerPitch = helikopterPart2_P.KalibrerPitch_Gain *
      rtb_HILReadEncoder_o2;
  }

  /* Integrator: '<S1>/Integrator'
   *
   * Regarding '<S1>/Integrator':
   *  Limited Integrator
   */
  if (helikopterPart2_X.Integrator_CSTATE >=
      helikopterPart2_P.Integrator_UpperSat ) {
    helikopterPart2_X.Integrator_CSTATE = helikopterPart2_P.Integrator_UpperSat;
  } else if (helikopterPart2_X.Integrator_CSTATE <=
             helikopterPart2_P.Integrator_LowerSat ) {
    helikopterPart2_X.Integrator_CSTATE = helikopterPart2_P.Integrator_LowerSat;
  }

  rtb_Saturation_h = helikopterPart2_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    /* Gain: '<S2>/Kalibrer -Vandring' */
    helikopterPart2_B.KalibrerVandring = helikopterPart2_P.KalibrerVandring_Gain
      * rtb_HILReadEncoder_o1;
  }

  /* TransferFcn: '<S2>/Vandring Deriv' */
  rtb_VandringDeriv = helikopterPart2_P.VandringDeriv_D*
    helikopterPart2_B.KalibrerVandring;
  rtb_VandringDeriv += helikopterPart2_P.VandringDeriv_C*
    helikopterPart2_X.VandringDeriv_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn4' */
  rtb_Gain2 = helikopterPart2_P.TransferFcn4_D*helikopterPart2_B.KalibrerPitch;
  rtb_Gain2 += helikopterPart2_P.TransferFcn4_C*
    helikopterPart2_X.TransferFcn4_CSTATE;

  /* TransferFcn: '<S2>/Transfer Fcn5' */
  rtb_Gain1 = helikopterPart2_P.TransferFcn5_D*helikopterPart2_B.KalibrerElev;
  rtb_Gain1 += helikopterPart2_P.TransferFcn5_C*
    helikopterPart2_X.TransferFcn5_CSTATE;

  /* Gain: '<Root>/Gain' incorporates:
   *  SignalConversion: '<Root>/TmpSignal ConversionAtGainInport1'
   */
  tmp[0] = helikopterPart2_B.VandringLavpass;
  tmp[1] = rtb_VandringDeriv;
  tmp[2] = helikopterPart2_B.KalibrerPitch;
  tmp[3] = rtb_Gain2;
  tmp[4] = helikopterPart2_B.Add;
  tmp[5] = rtb_Gain1;
  for (tmp_0 = 0; tmp_0 < 6; tmp_0++) {
    rtb_Gain[tmp_0] = 0.0;
    for (tmp_1 = 0; tmp_1 < 6; tmp_1++) {
      rtb_Gain[tmp_0] += helikopterPart2_P.Gain_Gain[6 * tmp_1 + tmp_0] *
        tmp[tmp_1];
    }
  }

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<Root>/elevation'
   */
  rtb_Gain1 = helikopterPart2_P.elevation_Value - rtb_Gain[4];

  /* Gain: '<S1>/K_ei' */
  helikopterPart2_B.K_ei = helikopterPart2_P.K_ei_Gain * rtb_Gain1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Gain: '<S1>/K_ed'
   *  Gain: '<S1>/K_ep'
   */
  rtb_Saturation_h = (helikopterPart2_P.K_ep_Gain * rtb_Gain1 + rtb_Saturation_h)
    + helikopterPart2_P.K_ed_Gain * rtb_Gain[5];

  /* Saturate: '<S1>/Saturation' */
  rtb_Saturation_h = rt_SATURATE(rtb_Saturation_h,
    helikopterPart2_P.Saturation_LowerSat, helikopterPart2_P.Saturation_UpperSat);

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helikopterPart2_DWork.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopterPart2_DWork.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopterPart2_DWork.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopterPart2_M->Timing.t[0];

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

    helikopterPart2_DWork.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

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

  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/K_pd'
   *  Gain: '<S3>/K_pp'
   *  Saturate: '<S3>/Saturation'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Gain1 = (rt_SATURATE(rtb_Gain1, helikopterPart2_P.Saturation_LowerSat_h,
    helikopterPart2_P.Saturation_UpperSat_o) - rtb_Gain[2]) *
    helikopterPart2_P.K_pp_Gain - helikopterPart2_P.K_pd_Gain * rtb_Gain[3];

  /* Gain: '<S4>/Gain2' incorporates:
   *  Sum: '<S4>/Sum4'
   */
  rtb_Gain2 = (rtb_Saturation_h - rtb_Gain1) * helikopterPart2_P.Gain2_Gain;

  /* Saturate: '<S2>/Sat B' */
  helikopterPart2_B.SatB = rt_SATURATE(rtb_Gain2,
    helikopterPart2_P.SatB_LowerSat, helikopterPart2_P.SatB_UpperSat);
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
  }

  /* Gain: '<S4>/Gain1' incorporates:
   *  Sum: '<S4>/Sum3'
   */
  rtb_Gain1 = (rtb_Gain1 + rtb_Saturation_h) * helikopterPart2_P.Gain1_Gain;

  /* Saturate: '<S2>/Sat' */
  helikopterPart2_B.Sat = rt_SATURATE(rtb_Gain1, helikopterPart2_P.Sat_LowerSat,
    helikopterPart2_P.Sat_UpperSat);
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    /* S-Function (hil_write_analog_block): '<S2>/HIL Write Analog' */

    /* S-Function Block: helikopterPart2/Heli 3D/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopterPart2_DWork.HILWriteAnalog_Buffer[0] = helikopterPart2_B.SatB;
      helikopterPart2_DWork.HILWriteAnalog_Buffer[1] = helikopterPart2_B.Sat;
      result = hil_write_analog(helikopterPart2_DWork.HILInitialize_Card,
        helikopterPart2_P.HILWriteAnalog_Channels, 2,
        &helikopterPart2_DWork.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
      }
    }
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Model update function */
void helikopterPart2_update(int_T tid)
{
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    rt_ertODEUpdateContinuousStates(&helikopterPart2_M->solverInfo);
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
  if (!(++helikopterPart2_M->Timing.clockTick0)) {
    ++helikopterPart2_M->Timing.clockTickH0;
  }

  helikopterPart2_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopterPart2_M->solverInfo);
  if (rtmIsMajorTimeStep(helikopterPart2_M)) {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopterPart2_M->Timing.clockTick1)) {
      ++helikopterPart2_M->Timing.clockTickH1;
    }

    helikopterPart2_M->Timing.t[1] = helikopterPart2_M->Timing.clockTick1 *
      helikopterPart2_M->Timing.stepSize1 +
      helikopterPart2_M->Timing.clockTickH1 *
      helikopterPart2_M->Timing.stepSize1 * 4294967296.0;
  }

  /* tid is required for a uniform function interface.
   * Argument tid is not used in the function. */
  UNUSED_PARAMETER(tid);
}

/* Derivatives for root system: '<Root>' */
void helikopterPart2_derivatives(void)
{
  /* Derivatives for TransferFcn: '<S2>/Vandring Lavpass' */
  {
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->VandringLavpass_CSTATE = helikopterPart2_B.KalibrerVandring;
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->VandringLavpass_CSTATE += (helikopterPart2_P.VandringLavpass_A)*
      helikopterPart2_X.VandringLavpass_CSTATE;
  }

  /* Derivatives for Integrator: '<S1>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopterPart2_X.Integrator_CSTATE <=
            helikopterPart2_P.Integrator_LowerSat );
    usat = ( helikopterPart2_X.Integrator_CSTATE >=
            helikopterPart2_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopterPart2_B.K_ei > 0)) ||
        (usat && (helikopterPart2_B.K_ei < 0)) ) {
      ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs
        )->Integrator_CSTATE = helikopterPart2_B.K_ei;
    } else {
      /* in saturation */
      ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs
        )->Integrator_CSTATE = 0.0;
    }
  }

  /* Derivatives for TransferFcn: '<S2>/Vandring Deriv' */
  {
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->VandringDeriv_CSTATE = helikopterPart2_B.KalibrerVandring;
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->VandringDeriv_CSTATE += (helikopterPart2_P.VandringDeriv_A)*
      helikopterPart2_X.VandringDeriv_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn4' */
  {
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->TransferFcn4_CSTATE = helikopterPart2_B.KalibrerPitch;
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->TransferFcn4_CSTATE += (helikopterPart2_P.TransferFcn4_A)*
      helikopterPart2_X.TransferFcn4_CSTATE;
  }

  /* Derivatives for TransferFcn: '<S2>/Transfer Fcn5' */
  {
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->TransferFcn5_CSTATE = helikopterPart2_B.KalibrerElev;
    ((StateDerivatives_helikopterPart *) helikopterPart2_M->ModelData.derivs)
      ->TransferFcn5_CSTATE += (helikopterPart2_P.TransferFcn5_A)*
      helikopterPart2_X.TransferFcn5_CSTATE;
  }
}

/* Model initialize function */
void helikopterPart2_initialize(boolean_T firstTime)
{
  (void)firstTime;

  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopterPart2_P.Integrator_UpperSat = rtInf;
  helikopterPart2_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopterPart2_M, 0,
                sizeof(RT_MODEL_helikopterPart2));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopterPart2_M->solverInfo,
                          &helikopterPart2_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopterPart2_M->solverInfo, &rtmGetTPtr(helikopterPart2_M));
    rtsiSetStepSizePtr(&helikopterPart2_M->solverInfo,
                       &helikopterPart2_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopterPart2_M->solverInfo,
                 &helikopterPart2_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopterPart2_M->solverInfo,
                         &helikopterPart2_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopterPart2_M->solverInfo,
      &helikopterPart2_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopterPart2_M->solverInfo, (&rtmGetErrorStatus
      (helikopterPart2_M)));
    rtsiSetRTModelPtr(&helikopterPart2_M->solverInfo, helikopterPart2_M);
  }

  rtsiSetSimTimeStep(&helikopterPart2_M->solverInfo, MAJOR_TIME_STEP);
  helikopterPart2_M->ModelData.intgData.f[0] = helikopterPart2_M->
    ModelData.odeF[0];
  helikopterPart2_M->ModelData.contStates = ((real_T *) &helikopterPart2_X);
  rtsiSetSolverData(&helikopterPart2_M->solverInfo, (void *)
                    &helikopterPart2_M->ModelData.intgData);
  rtsiSetSolverName(&helikopterPart2_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopterPart2_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopterPart2_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopterPart2_M->Timing.sampleTimes =
      (&helikopterPart2_M->Timing.sampleTimesArray[0]);
    helikopterPart2_M->Timing.offsetTimes =
      (&helikopterPart2_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopterPart2_M->Timing.sampleTimes[0] = (0.0);
    helikopterPart2_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    helikopterPart2_M->Timing.offsetTimes[0] = (0.0);
    helikopterPart2_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopterPart2_M, &helikopterPart2_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopterPart2_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopterPart2_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopterPart2_M, -1);
  helikopterPart2_M->Timing.stepSize0 = 0.001;
  helikopterPart2_M->Timing.stepSize1 = 0.001;

  /* external mode info */
  helikopterPart2_M->Sizes.checksums[0] = (962577980U);
  helikopterPart2_M->Sizes.checksums[1] = (1449715909U);
  helikopterPart2_M->Sizes.checksums[2] = (1724484102U);
  helikopterPart2_M->Sizes.checksums[3] = (927469709U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopterPart2_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopterPart2_M->extModeInfo,
      &helikopterPart2_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopterPart2_M->extModeInfo,
                        helikopterPart2_M->Sizes.checksums);
    rteiSetTPtr(helikopterPart2_M->extModeInfo, rtmGetTPtr(helikopterPart2_M));
  }

  helikopterPart2_M->solverInfoPtr = (&helikopterPart2_M->solverInfo);
  helikopterPart2_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&helikopterPart2_M->solverInfo, 0.001);
  rtsiSetSolverMode(&helikopterPart2_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopterPart2_M->ModelData.blockIO = ((void *) &helikopterPart2_B);

  {
    helikopterPart2_B.KalibrerElev = 0.0;
    helikopterPart2_B.Add = 0.0;
    helikopterPart2_B.VandringLavpass = 0.0;
    helikopterPart2_B.KalibrerPitch = 0.0;
    helikopterPart2_B.KalibrerVandring = 0.0;
    helikopterPart2_B.K_ei = 0.0;
    helikopterPart2_B.SatB = 0.0;
    helikopterPart2_B.Sat = 0.0;
  }

  /* parameters */
  helikopterPart2_M->ModelData.defaultParam = ((real_T *)&helikopterPart2_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopterPart2_X;
    helikopterPart2_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopterPart2_X, 0,
                  sizeof(ContinuousStates_helikopterPart));
  }

  /* states (dwork) */
  helikopterPart2_M->Work.dwork = ((void *) &helikopterPart2_DWork);
  (void) memset((void *)&helikopterPart2_DWork, 0,
                sizeof(D_Work_helikopterPart2));

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterPart2_DWork.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 16; i++) {
      helikopterPart2_DWork.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  helikopterPart2_DWork.HILInitialize_AOMinimums[0] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMinimums[1] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMinimums[2] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMinimums[3] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMaximums[0] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMaximums[1] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMaximums[2] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOMaximums[3] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOVoltages[0] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOVoltages[1] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOVoltages[2] = 0.0;
  helikopterPart2_DWork.HILInitialize_AOVoltages[3] = 0.0;
  helikopterPart2_DWork.HILWriteAnalog_Buffer[0] = 0.0;
  helikopterPart2_DWork.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopterPart2_M->SpecialInfo.mappingInfo = (&dtInfo);
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
void helikopterPart2_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterPart2/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    hil_task_stop_all(helikopterPart2_DWork.HILInitialize_Card);
    hil_task_delete_all(helikopterPart2_DWork.HILInitialize_Card);
    hil_monitor_stop_all(helikopterPart2_DWork.HILInitialize_Card);
    hil_monitor_delete_all(helikopterPart2_DWork.HILInitialize_Card);
    is_switching = false;
    if ((helikopterPart2_P.HILInitialize_AOTerminate && !is_switching) ||
        (helikopterPart2_P.HILInitialize_AOExit && is_switching)) {
      helikopterPart2_DWork.HILInitialize_AOVoltages[0] =
        helikopterPart2_P.HILInitialize_AOFinal;
      helikopterPart2_DWork.HILInitialize_AOVoltages[1] =
        helikopterPart2_P.HILInitialize_AOFinal;
      helikopterPart2_DWork.HILInitialize_AOVoltages[2] =
        helikopterPart2_P.HILInitialize_AOFinal;
      helikopterPart2_DWork.HILInitialize_AOVoltages[3] =
        helikopterPart2_P.HILInitialize_AOFinal;
      result = hil_write_analog(helikopterPart2_DWork.HILInitialize_Card,
        &helikopterPart2_P.HILInitialize_AOChannels[0], 4U,
        &helikopterPart2_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
      }
    }

    hil_close(helikopterPart2_DWork.HILInitialize_Card);
    helikopterPart2_DWork.HILInitialize_Card = NULL;
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
  helikopterPart2_output(tid);
}

void MdlUpdate(int_T tid)
{
  helikopterPart2_update(tid);
}

void MdlInitializeSizes(void)
{
  helikopterPart2_M->Sizes.numContStates = (5);/* Number of continuous states */
  helikopterPart2_M->Sizes.numY = (0); /* Number of model outputs */
  helikopterPart2_M->Sizes.numU = (0); /* Number of model inputs */
  helikopterPart2_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopterPart2_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopterPart2_M->Sizes.numBlocks = (45);/* Number of blocks */
  helikopterPart2_M->Sizes.numBlockIO = (8);/* Number of block outputs */
  helikopterPart2_M->Sizes.numBlockPrms = (149);/* Sum of parameter "widths" */
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
  /* InitializeConditions for TransferFcn: '<S2>/Vandring Lavpass' */
  helikopterPart2_X.VandringLavpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S1>/Integrator' */
  helikopterPart2_X.Integrator_CSTATE = helikopterPart2_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S2>/Vandring Deriv' */
  helikopterPart2_X.VandringDeriv_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn4' */
  helikopterPart2_X.TransferFcn4_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S2>/Transfer Fcn5' */
  helikopterPart2_X.TransferFcn5_CSTATE = 0.0;
}

void MdlStart(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopterPart2/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("sensoray_model_626", "0",
                      &helikopterPart2_DWork.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
      return;
    }

    is_switching = false;
    if ((helikopterPart2_P.HILInitialize_CKPStart && !is_switching) ||
        (helikopterPart2_P.HILInitialize_CKPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_ClockModes =
          &helikopterPart2_DWork.HILInitialize_ClockModes[0];
        for (i1=0; i1 < 6; i1++) {
          dw_ClockModes[i1] = helikopterPart2_P.HILInitialize_CKModes;
        }
      }

      result = hil_set_clock_mode(helikopterPart2_DWork.HILInitialize_Card,
        (t_clock *) &helikopterPart2_P.HILInitialize_CKChannels[0], 6U,
        (t_clock_mode *) &helikopterPart2_DWork.HILInitialize_ClockModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
        return;
      }
    }

    result = hil_watchdog_clear(helikopterPart2_DWork.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
      return;
    }

    if ((helikopterPart2_P.HILInitialize_AIPStart && !is_switching) ||
        (helikopterPart2_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helikopterPart2_DWork.HILInitialize_AIMinimums
          [0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMinimums[i1] = helikopterPart2_P.HILInitialize_AILow;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helikopterPart2_DWork.HILInitialize_AIMaximums
          [0];
        for (i1=0; i1 < 16; i1++) {
          dw_AIMaximums[i1] = helikopterPart2_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (helikopterPart2_DWork.HILInitialize_Card,
         &helikopterPart2_P.HILInitialize_AIChannels[0], 16U,
         &helikopterPart2_DWork.HILInitialize_AIMinimums[0],
         &helikopterPart2_DWork.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterPart2_P.HILInitialize_AOPStart && !is_switching) ||
        (helikopterPart2_P.HILInitialize_AOPEnter && is_switching)) {
      helikopterPart2_DWork.HILInitialize_AOMinimums[0] =
        helikopterPart2_P.HILInitialize_AOLow;
      helikopterPart2_DWork.HILInitialize_AOMinimums[1] =
        helikopterPart2_P.HILInitialize_AOLow;
      helikopterPart2_DWork.HILInitialize_AOMinimums[2] =
        helikopterPart2_P.HILInitialize_AOLow;
      helikopterPart2_DWork.HILInitialize_AOMinimums[3] =
        helikopterPart2_P.HILInitialize_AOLow;
      helikopterPart2_DWork.HILInitialize_AOMaximums[0] =
        helikopterPart2_P.HILInitialize_AOHigh;
      helikopterPart2_DWork.HILInitialize_AOMaximums[1] =
        helikopterPart2_P.HILInitialize_AOHigh;
      helikopterPart2_DWork.HILInitialize_AOMaximums[2] =
        helikopterPart2_P.HILInitialize_AOHigh;
      helikopterPart2_DWork.HILInitialize_AOMaximums[3] =
        helikopterPart2_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (helikopterPart2_DWork.HILInitialize_Card,
         &helikopterPart2_P.HILInitialize_AOChannels[0], 4U,
         &helikopterPart2_DWork.HILInitialize_AOMinimums[0],
         &helikopterPart2_DWork.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
        return;
      }
    }

    if ((helikopterPart2_P.HILInitialize_AOStart && !is_switching) ||
        (helikopterPart2_P.HILInitialize_AOEnter && is_switching)) {
      helikopterPart2_DWork.HILInitialize_AOVoltages[0] =
        helikopterPart2_P.HILInitialize_AOInitial;
      helikopterPart2_DWork.HILInitialize_AOVoltages[1] =
        helikopterPart2_P.HILInitialize_AOInitial;
      helikopterPart2_DWork.HILInitialize_AOVoltages[2] =
        helikopterPart2_P.HILInitialize_AOInitial;
      helikopterPart2_DWork.HILInitialize_AOVoltages[3] =
        helikopterPart2_P.HILInitialize_AOInitial;
      result = hil_write_analog(helikopterPart2_DWork.HILInitialize_Card,
        &helikopterPart2_P.HILInitialize_AOChannels[0], 4U,
        &helikopterPart2_DWork.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopterPart2_M, _rt_error_message);
        return;
      }
    }
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
      1.6455121993178662E-001, 8.2275609965890867E-002, 5.2359877559829882E-001,
      5.2359877559829460E-001, 5.2359877559829782E-001, 5.2359877559829715E-001,
      5.2359877559829771E-001, 5.2359877559829782E-001, 5.2359877559829882E-001,
      5.2359877559829882E-001, 5.2359877559829882E-001, 5.2359877559829759E-001,
      5.2359877559829882E-001, 5.2359877559829859E-001, 5.2359877559829870E-001,
      5.2359877559829882E-001, 5.2359877559829449E-001, 5.2359877559829582E-001,
      5.2359877559829548E-001, 5.2359877559829604E-001, 5.2359877559829859E-001,
      5.2359877559829759E-001, 5.2359877559829837E-001, 5.2359877559829504E-001,
      5.2359877559829882E-001, 5.1054365978924121E-001, 4.6255718828140246E-001,
      4.0234351749718478E-001, 2.5042947270242671E-001, 1.0126211051142006E-001,
      1.3145866937202583E-003, -7.6739201443561292E-002,
      -1.5609842989698350E-001, -2.3167410763544863E-001,
      -2.9463551283375100E-001, -3.4599057507076036E-001,
      -3.8996159150599446E-001, -4.2768369780644000E-001,
      -4.5849738609861984E-001, -4.8253676826668918E-001,
      -5.0072831920997340E-001, -5.1382414193405490E-001,
      -5.2222105402904084E-001, -5.2359877559829737E-001,
      -5.2333250385733820E-001, -5.2104484084821900E-001,
      -5.1784170054901080E-001, -5.0921828028514704E-001,
      -4.9700566078136849E-001, -4.8362569892323110E-001,
      -4.6922019672515880E-001, -4.5325850748665569E-001,
      -4.3590840797768660E-001, -4.1769311085749794E-001,
      -3.9891711968409987E-001, -3.7968806898999374E-001,
      -3.6013589634141091E-001, -3.4044661432117773E-001,
      -3.2078796434361512E-001, -3.0128405348630066E-001,
      -2.8203943728672770E-001, -2.6315415121296548E-001,
      -2.4471860663337056E-001, -2.2680809019990880E-001,
      -2.0948467399994522E-001, -1.9280062772060425E-001,
      -1.7679931654673045E-001, -1.6151503604535750E-001,
      -1.4697354946041311E-001, -1.3319311332042147E-001,
      -1.2018529157628242E-001, -1.0795550647239552E-001,
      -9.6503576629527313E-002, -8.5824307875172298E-002,
      -7.5908055478616954E-002, -6.6741218428576901E-002,
      -5.8306689339102373E-002, -5.0584279635219200E-002,
      -4.3551115271599668E-002, -3.7181995810101870E-002,
      -3.1449719780833246E-002, -2.6325381562437404E-002,
      -2.1778638007289533E-002, -1.7777940727723499E-002,
      -1.4290742270999558E-002, -1.1283693321250697E-002,
      -8.7228210885951420E-003, -6.5736434643052424E-003,
      -4.8012279903062016E-003, -3.3703268186026526E-003,
      -2.2456338459954372E-003, -1.3918404790898649E-003,
      -7.7320967359090686E-004, -3.5340107010761935E-004,
      -9.6670757346880624E-005, 3.0952764430801324E-005, 6.4223612444222464E-005,
      4.3062150628388074E-005, 7.8680910188487079E-006, 7.6586349095262388E-015,
      7.6586349095262388E-015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopterPart2_DWork.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues;
    helikopterPart2_DWork.FromWorkspace_PWORK.DataPtr = (void *) pDataValues;
    helikopterPart2_DWork.FromWorkspace_IWORK.PrevIndex = 0;
  }

  MdlInitialize();
}

void MdlTerminate(void)
{
  helikopterPart2_terminate();
}

RT_MODEL_helikopterPart2 *helikopterPart2(void)
{
  helikopterPart2_initialize(1);
  return helikopterPart2_M;
}

/*========================================================================*
 * End of GRT compatible call interface                                   *
 *========================================================================*/
