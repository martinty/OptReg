/*
 * helicopter2.c
 *
 * Code generation for model "helicopter2".
 *
 * Model version              : 1.179
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Apr 27 12:42:45 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter2.h"
#include "helicopter2_private.h"
#include "helicopter2_dt.h"

/* Block signals (auto storage) */
B_helicopter2_T helicopter2_B;

/* Continuous states */
X_helicopter2_T helicopter2_X;

/* Block states (auto storage) */
DW_helicopter2_T helicopter2_DW;

/* Real-time model */
RT_MODEL_helicopter2_T helicopter2_M_;
RT_MODEL_helicopter2_T *const helicopter2_M = &helicopter2_M_;

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
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter2_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter2_output(void)
{
  /* local block i/o variables */
  real_T rtb_Frontgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Backgain;
  real_T rtb_Derivative;
  real_T rtb_Gain1_idx_4;
  real_T rtb_Gain1_idx_5;
  if (rtmIsMajorTimeStep(helicopter2_M)) {
    /* set solver stop time */
    if (!(helicopter2_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter2_M->solverInfo,
                            ((helicopter2_M->Timing.clockTickH0 + 1) *
        helicopter2_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter2_M->solverInfo,
                            ((helicopter2_M->Timing.clockTick0 + 1) *
        helicopter2_M->Timing.stepSize0 + helicopter2_M->Timing.clockTickH0 *
        helicopter2_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter2_M)) {
    helicopter2_M->Timing.t[0] = rtsiGetT(&helicopter2_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter2_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S6>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter2/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter2_DW.HILReadEncoderTimebase_Task,
        1, &helicopter2_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter2_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter2_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter2_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S6>/Travel: Count to rad' */
    helicopter2_B.TravelCounttorad = helicopter2_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S13>/Gain' */
    helicopter2_B.Gain = helicopter2_P.Gain_Gain *
      helicopter2_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helicopter2_B.Sum1 = helicopter2_P.Constant_Value + helicopter2_B.Gain;

    /* Gain: '<S6>/Pitch: Count to rad' */
    helicopter2_B.PitchCounttorad = helicopter2_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S10>/Gain' */
    helicopter2_B.Gain_i = helicopter2_P.Gain_Gain_a *
      helicopter2_B.PitchCounttorad;

    /* Gain: '<S4>/Gain1' */
    helicopter2_B.Gain1[0] = helicopter2_P.Gain1_Gain * helicopter2_B.Sum1;
    helicopter2_B.Gain1[1] = helicopter2_P.Gain1_Gain * helicopter2_B.Gain_i;
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *) helicopter2_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter2_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter2_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter2_M->Timing.t[0];

    /* Get index */
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

    helicopter2_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Frontgain = pDataValues[currTimeIndex];
        } else {
          rtb_Frontgain = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Frontgain = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter2_M)) {
    /* Gain: '<S3>/Gain1' incorporates:
     *  Constant: '<Root>/Constant1'
     */
    helicopter2_B.Gain1_a = helicopter2_P.Gain1_Gain_a *
      helicopter2_P.Constant1_Value;
  }

  /* Gain: '<S14>/Gain' incorporates:
   *  TransferFcn: '<S6>/Travel: Transfer Fcn'
   */
  helicopter2_B.Gain_d = (helicopter2_P.TravelTransferFcn_C *
    helicopter2_X.TravelTransferFcn_CSTATE + helicopter2_P.TravelTransferFcn_D *
    helicopter2_B.TravelCounttorad) * helicopter2_P.Gain_Gain_l;

  /* Gain: '<S11>/Gain' incorporates:
   *  TransferFcn: '<S6>/Pitch: Transfer Fcn'
   */
  helicopter2_B.Gain_b = (helicopter2_P.PitchTransferFcn_C *
    helicopter2_X.PitchTransferFcn_CSTATE + helicopter2_P.PitchTransferFcn_D *
    helicopter2_B.PitchCounttorad) * helicopter2_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(helicopter2_M)) {
    /* Gain: '<S6>/Elevation: Count to rad' */
    helicopter2_B.ElevationCounttorad = helicopter2_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S8>/Gain' */
    helicopter2_B.Gain_e = helicopter2_P.Gain_Gain_lv *
      helicopter2_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter2_B.Sum = helicopter2_B.Gain_e +
      helicopter2_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S6>/Elevation: Transfer Fcn'
   */
  helicopter2_B.Gain_dg = (helicopter2_P.ElevationTransferFcn_C *
    helicopter2_X.ElevationTransferFcn_CSTATE +
    helicopter2_P.ElevationTransferFcn_D * helicopter2_B.ElevationCounttorad) *
    helicopter2_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helicopter2_P.Gain1_Gain_f * helicopter2_B.Sum;
  rtb_Gain1_idx_5 = helicopter2_P.Gain1_Gain_f * helicopter2_B.Gain_dg;

  /* Sum: '<S7>/Sum' incorporates:
   *  Constant: '<S7>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S7>/K_pd'
   *  Gain: '<S7>/K_pp'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<S7>/Sum2'
   *  Sum: '<S7>/Sum3'
   */
  rtb_Backgain = ((rtb_Frontgain - (helicopter2_P.Gain1_Gain_f *
    helicopter2_B.Gain_i + helicopter2_B.Gain1_a)) * helicopter2_P.K_pp -
                  helicopter2_P.Gain1_Gain_f * helicopter2_B.Gain_b *
                  helicopter2_P.K_pd) + helicopter2_P.Vd_ff;

  /* Integrator: '<S5>/Integrator'
   *
   * Regarding '<S5>/Integrator':
   *  Limited Integrator
   */
  if (helicopter2_X.Integrator_CSTATE >= helicopter2_P.Integrator_UpperSat ) {
    helicopter2_X.Integrator_CSTATE = helicopter2_P.Integrator_UpperSat;
  } else if (helicopter2_X.Integrator_CSTATE <=
             (helicopter2_P.Integrator_LowerSat) ) {
    helicopter2_X.Integrator_CSTATE = (helicopter2_P.Integrator_LowerSat);
  }

  rtb_Frontgain = helicopter2_X.Integrator_CSTATE;

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helicopter2_P.elevation_ref_Value - rtb_Gain1_idx_4;

  /* Sum: '<S5>/Sum2' incorporates:
   *  Constant: '<S5>/Vs_bias'
   *  Gain: '<S5>/K_ed'
   *  Gain: '<S5>/K_ep'
   *  Sum: '<S5>/Sum1'
   */
  rtb_Frontgain = ((helicopter2_P.K_ep * rtb_Derivative + rtb_Frontgain) -
                   helicopter2_P.K_ed * rtb_Gain1_idx_5) + helicopter2_P.Vs_ff;

  /* Sum: '<S1>/Subtract' */
  rtb_Gain1_idx_4 = rtb_Frontgain - rtb_Backgain;

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Frontgain = (rtb_Backgain + rtb_Frontgain) * helicopter2_P.Frontgain_Gain;

  /* Gain: '<S5>/K_ei' */
  helicopter2_B.K_ei = helicopter2_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter2_M)) {
  }

  /* Derivative: '<S6>/Derivative' */
  if ((helicopter2_DW.TimeStampA >= helicopter2_M->Timing.t[0]) &&
      (helicopter2_DW.TimeStampB >= helicopter2_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Backgain = helicopter2_DW.TimeStampA;
    lastU = &helicopter2_DW.LastUAtTimeA;
    if (helicopter2_DW.TimeStampA < helicopter2_DW.TimeStampB) {
      if (helicopter2_DW.TimeStampB < helicopter2_M->Timing.t[0]) {
        rtb_Backgain = helicopter2_DW.TimeStampB;
        lastU = &helicopter2_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter2_DW.TimeStampA >= helicopter2_M->Timing.t[0]) {
        rtb_Backgain = helicopter2_DW.TimeStampB;
        lastU = &helicopter2_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter2_B.PitchCounttorad - *lastU) /
      (helicopter2_M->Timing.t[0] - rtb_Backgain);
  }

  /* End of Derivative: '<S6>/Derivative' */

  /* Gain: '<S12>/Gain' */
  helicopter2_B.Gain_l = helicopter2_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter2_M)) {
  }

  /* Gain: '<S1>/Back gain' */
  rtb_Gain1_idx_4 *= helicopter2_P.Backgain_Gain;

  /* Saturate: '<S6>/Back motor: Saturation' */
  if (rtb_Gain1_idx_4 > helicopter2_P.BackmotorSaturation_UpperSat) {
    helicopter2_B.BackmotorSaturation =
      helicopter2_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 < helicopter2_P.BackmotorSaturation_LowerSat) {
    helicopter2_B.BackmotorSaturation =
      helicopter2_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter2_B.BackmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S6>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter2_M)) {
  }

  /* Saturate: '<S6>/Front motor: Saturation' */
  if (rtb_Frontgain > helicopter2_P.FrontmotorSaturation_UpperSat) {
    helicopter2_B.FrontmotorSaturation =
      helicopter2_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Frontgain < helicopter2_P.FrontmotorSaturation_LowerSat) {
    helicopter2_B.FrontmotorSaturation =
      helicopter2_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter2_B.FrontmotorSaturation = rtb_Frontgain;
  }

  /* End of Saturate: '<S6>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter2_M)) {
    /* S-Function (hil_write_analog_block): '<S6>/HIL Write Analog' */

    /* S-Function Block: helicopter2/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter2_DW.HILWriteAnalog_Buffer[0] =
        helicopter2_B.FrontmotorSaturation;
      helicopter2_DW.HILWriteAnalog_Buffer[1] =
        helicopter2_B.BackmotorSaturation;
      result = hil_write_analog(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILWriteAnalog_channels, 2,
        &helicopter2_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter2_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S6>/Derivative' */
  if (helicopter2_DW.TimeStampA == (rtInf)) {
    helicopter2_DW.TimeStampA = helicopter2_M->Timing.t[0];
    lastU = &helicopter2_DW.LastUAtTimeA;
  } else if (helicopter2_DW.TimeStampB == (rtInf)) {
    helicopter2_DW.TimeStampB = helicopter2_M->Timing.t[0];
    lastU = &helicopter2_DW.LastUAtTimeB;
  } else if (helicopter2_DW.TimeStampA < helicopter2_DW.TimeStampB) {
    helicopter2_DW.TimeStampA = helicopter2_M->Timing.t[0];
    lastU = &helicopter2_DW.LastUAtTimeA;
  } else {
    helicopter2_DW.TimeStampB = helicopter2_M->Timing.t[0];
    lastU = &helicopter2_DW.LastUAtTimeB;
  }

  *lastU = helicopter2_B.PitchCounttorad;

  /* End of Update for Derivative: '<S6>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter2_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter2_M->solverInfo);
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
  if (!(++helicopter2_M->Timing.clockTick0)) {
    ++helicopter2_M->Timing.clockTickH0;
  }

  helicopter2_M->Timing.t[0] = rtsiGetSolverStopTime(&helicopter2_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter2_M->Timing.clockTick1)) {
      ++helicopter2_M->Timing.clockTickH1;
    }

    helicopter2_M->Timing.t[1] = helicopter2_M->Timing.clockTick1 *
      helicopter2_M->Timing.stepSize1 + helicopter2_M->Timing.clockTickH1 *
      helicopter2_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter2_derivatives(void)
{
  XDot_helicopter2_T *_rtXdot;
  _rtXdot = ((XDot_helicopter2_T *) helicopter2_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S6>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter2_P.TravelTransferFcn_A *
    helicopter2_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter2_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S6>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter2_P.PitchTransferFcn_A *
    helicopter2_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter2_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S6>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter2_P.ElevationTransferFcn_A *
    helicopter2_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter2_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S5>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter2_X.Integrator_CSTATE <=
            (helicopter2_P.Integrator_LowerSat) );
    usat = ( helicopter2_X.Integrator_CSTATE >=
            helicopter2_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter2_B.K_ei > 0)) ||
        (usat && (helicopter2_B.K_ei < 0)) ) {
      ((XDot_helicopter2_T *) helicopter2_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter2_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter2_T *) helicopter2_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter2_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter2/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter2_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter2_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter2_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter2_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter2_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter2_M, _rt_error_message);
      return;
    }

    if ((helicopter2_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter2_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter2_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter2_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter2_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter2_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_analog_input_chan, 8U,
        &helicopter2_DW.HILInitialize_AIMinimums[0],
        &helicopter2_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter2_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter2_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter2_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter2_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter2_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter2_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_analog_output_cha, 8U,
        &helicopter2_DW.HILInitialize_AOMinimums[0],
        &helicopter2_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter2_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter2_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter2_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter2_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_analog_output_cha, 8U,
        &helicopter2_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if (helicopter2_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter2_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter2_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter2_DW.HILInitialize_Card,
         helicopter2_P.HILInitialize_analog_output_cha, 8U,
         &helicopter2_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter2_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter2_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter2_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter2_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_encoder_channels, 8U,
        (t_encoder_quadrature_mode *)
        &helicopter2_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter2_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter2_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter2_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicopter2_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_encoder_channels, 8U,
        &helicopter2_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter2_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter2_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter2_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter2_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter2_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter2_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helicopter2_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter2_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter2_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter2_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter2_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter2_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helicopter2_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter2_DW.HILInitialize_Card,
          &helicopter2_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter2_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter2_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter2_DW.HILInitialize_Card,
          &helicopter2_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter2_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter2_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter2_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter2_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &helicopter2_DW.HILInitialize_POAlignValues
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter2_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter2_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter2_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter2_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter2_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter2_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter2_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter2_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter2_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_pwm_channels, 8U,
        &helicopter2_DW.HILInitialize_POSortedFreqs[0],
        &helicopter2_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter2_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter2_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter2_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter2_DW.HILInitialize_Card,
        helicopter2_P.HILInitialize_pwm_channels, 8U,
        &helicopter2_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }

    if (helicopter2_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter2_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter2_DW.HILInitialize_Card,
         helicopter2_P.HILInitialize_pwm_channels, 8U,
         &helicopter2_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S6>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter2/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter2_DW.HILInitialize_Card,
      helicopter2_P.HILReadEncoderTimebase_samples_,
      helicopter2_P.HILReadEncoderTimebase_channels, 3,
      &helicopter2_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter2_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
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

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559807444,
      0.52359877559803447, 0.52359877559798307, 0.52359877559791523,
      0.523598775597823, 0.5235987755976923, 0.52359877559749723,
      0.52359877559718215, 0.5235987755966055, 0.52359877559526735,
      0.52359877558925361, 0.38860546189573791, 0.10951698887526917,
      -0.1100383394610223, -0.27691123755291591, -0.3979068273164692,
      -0.47963675127901134, -0.523591909441963, -0.52359877534342569,
      -0.52359877536343258, -0.5235987224145312, -0.50343489391626584,
      -0.46497044043306424, -0.42077345861352544, -0.37347862555446432,
      -0.32522986550994049, -0.27772344682220712, -0.2322553027667397,
      -0.1897701881732447, -0.15091074088283754, -0.11606494657763659,
      -0.08541089068498342, -0.058958016011036223, -0.036584388943327412,
      -0.018069712651425172, -0.0031240162861197467, 0.008587900985857537,
      0.017426078578606438, 0.023758780967378875, 0.027949332008985642,
      0.030346034333594667, 0.03127479682777045, 0.0310340938210644,
      0.0298918916543043, 0.02808419891071846, 0.025814923189764548,
      0.023256747740391905, 0.020552773751695553, 0.017818707168261894,
      0.015145401405518497, 0.01260159841128787, 0.010236739517233149,
      0.0080837440173715247, 0.0061616771427812711, 0.0044782499569109863,
      0.0030321116735972076, 0.0018149100878760917, 0.00081310836210016267,
      9.5565276842160941E-6, -0.00061517602288165755, -0.0010816948208518151,
      -0.0014108848952404636, -0.0016232063484833872, -0.0017381600630382574,
      -0.0017739020793310287, -0.0017469853873323285, -0.0016722086841820883,
      -0.0015625529081204723, -0.0014291879259837978, -0.0012815335105175723,
      -0.0011273605978063502, -0.00097292068600022858, -0.00082309306297666392,
      -0.00068154128644839379, -0.00055087195231493555, -0.00043279025403042593,
      -0.00032824814515925653, -0.00023758206462032373, -0.00016063817112712349,
      -9.688386664573433E-5, -4.5505078677972685E-5, -5.4893309923764594E-6,
      2.4304922656237162E-5, 4.5091880270843108E-5, 5.8113799538615989E-5,
      6.4606933573472039E-5, 6.5776522273753571E-5, 6.2779201128075149E-5,
      5.6710998483779331E-5, 4.8598837981454231E-5, 3.9393103171466183E-5,
      2.9958323012006171E-5, 2.1058409934228412E-5, 1.3332266149043419E-5,
      7.2554391803702282E-6, 3.0850995699961848E-6, 7.9190498794152739E-7,
      -4.0593744014257015E-17, -1.5332144970613159E-17, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicopter2_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter2_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter2_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S6>/Travel: Transfer Fcn' */
  helicopter2_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S6>/Pitch: Transfer Fcn' */
  helicopter2_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S6>/Elevation: Transfer Fcn' */
  helicopter2_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S5>/Integrator' */
  helicopter2_X.Integrator_CSTATE = helicopter2_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S6>/Derivative' */
  helicopter2_DW.TimeStampA = (rtInf);
  helicopter2_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter2_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter2/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter2_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter2_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter2_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter2_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter2_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter2_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter2_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter2_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter2_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter2_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter2_DW.HILInitialize_Card
                         , helicopter2_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter2_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter2_DW.HILInitialize_AOVoltages[0]
                         , &helicopter2_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter2_DW.HILInitialize_Card,
            helicopter2_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helicopter2_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter2_DW.HILInitialize_Card,
            helicopter2_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter2_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter2_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter2_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter2_DW.HILInitialize_Card);
    hil_close(helicopter2_DW.HILInitialize_Card);
    helicopter2_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
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
  helicopter2_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter2_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter2_initialize();
}

void MdlTerminate(void)
{
  helicopter2_terminate();
}

/* Registration function */
RT_MODEL_helicopter2_T *helicopter2(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter2_P.Integrator_UpperSat = rtInf;
  helicopter2_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter2_M, 0,
                sizeof(RT_MODEL_helicopter2_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter2_M->solverInfo,
                          &helicopter2_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter2_M->solverInfo, &rtmGetTPtr(helicopter2_M));
    rtsiSetStepSizePtr(&helicopter2_M->solverInfo,
                       &helicopter2_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter2_M->solverInfo, &helicopter2_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter2_M->solverInfo, (real_T **)
                         &helicopter2_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter2_M->solverInfo,
      &helicopter2_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter2_M->solverInfo, (&rtmGetErrorStatus
      (helicopter2_M)));
    rtsiSetRTModelPtr(&helicopter2_M->solverInfo, helicopter2_M);
  }

  rtsiSetSimTimeStep(&helicopter2_M->solverInfo, MAJOR_TIME_STEP);
  helicopter2_M->ModelData.intgData.f[0] = helicopter2_M->ModelData.odeF[0];
  helicopter2_M->ModelData.contStates = ((real_T *) &helicopter2_X);
  rtsiSetSolverData(&helicopter2_M->solverInfo, (void *)
                    &helicopter2_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter2_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter2_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter2_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter2_M->Timing.sampleTimes = (&helicopter2_M->
      Timing.sampleTimesArray[0]);
    helicopter2_M->Timing.offsetTimes = (&helicopter2_M->
      Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter2_M->Timing.sampleTimes[0] = (0.0);
    helicopter2_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter2_M->Timing.offsetTimes[0] = (0.0);
    helicopter2_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter2_M, &helicopter2_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter2_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter2_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter2_M, 20.0);
  helicopter2_M->Timing.stepSize0 = 0.002;
  helicopter2_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter2_M->Sizes.checksums[0] = (3708942728U);
  helicopter2_M->Sizes.checksums[1] = (108629731U);
  helicopter2_M->Sizes.checksums[2] = (2953796357U);
  helicopter2_M->Sizes.checksums[3] = (2469943168U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter2_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter2_M->extModeInfo,
      &helicopter2_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter2_M->extModeInfo,
                        helicopter2_M->Sizes.checksums);
    rteiSetTPtr(helicopter2_M->extModeInfo, rtmGetTPtr(helicopter2_M));
  }

  helicopter2_M->solverInfoPtr = (&helicopter2_M->solverInfo);
  helicopter2_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter2_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter2_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter2_M->ModelData.blockIO = ((void *) &helicopter2_B);

  {
    helicopter2_B.TravelCounttorad = 0.0;
    helicopter2_B.Gain = 0.0;
    helicopter2_B.Sum1 = 0.0;
    helicopter2_B.PitchCounttorad = 0.0;
    helicopter2_B.Gain_i = 0.0;
    helicopter2_B.Gain1[0] = 0.0;
    helicopter2_B.Gain1[1] = 0.0;
    helicopter2_B.Gain1_a = 0.0;
    helicopter2_B.Gain_d = 0.0;
    helicopter2_B.Gain_b = 0.0;
    helicopter2_B.ElevationCounttorad = 0.0;
    helicopter2_B.Gain_e = 0.0;
    helicopter2_B.Sum = 0.0;
    helicopter2_B.Gain_dg = 0.0;
    helicopter2_B.K_ei = 0.0;
    helicopter2_B.Gain_l = 0.0;
    helicopter2_B.BackmotorSaturation = 0.0;
    helicopter2_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter2_M->ModelData.defaultParam = ((real_T *)&helicopter2_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter2_X;
    helicopter2_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter2_X, 0,
                  sizeof(X_helicopter2_T));
  }

  /* states (dwork) */
  helicopter2_M->ModelData.dwork = ((void *) &helicopter2_DW);
  (void) memset((void *)&helicopter2_DW, 0,
                sizeof(DW_helicopter2_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter2_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter2_DW.TimeStampA = 0.0;
  helicopter2_DW.LastUAtTimeA = 0.0;
  helicopter2_DW.TimeStampB = 0.0;
  helicopter2_DW.LastUAtTimeB = 0.0;
  helicopter2_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter2_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter2_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter2_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter2_M->Sizes.numY = (0);     /* Number of model outputs */
  helicopter2_M->Sizes.numU = (0);     /* Number of model inputs */
  helicopter2_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter2_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter2_M->Sizes.numBlocks = (58);/* Number of blocks */
  helicopter2_M->Sizes.numBlockIO = (17);/* Number of block outputs */
  helicopter2_M->Sizes.numBlockPrms = (145);/* Sum of parameter "widths" */
  return helicopter2_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
