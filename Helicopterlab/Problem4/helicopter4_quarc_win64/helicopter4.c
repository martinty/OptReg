/*
 * helicopter4.c
 *
 * Code generation for model "helicopter4".
 *
 * Model version              : 1.183
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Apr 27 15:16:35 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter4.h"
#include "helicopter4_private.h"
#include "helicopter4_dt.h"

/* Block signals (auto storage) */
B_helicopter4_T helicopter4_B;

/* Continuous states */
X_helicopter4_T helicopter4_X;

/* Block states (auto storage) */
DW_helicopter4_T helicopter4_DW;

/* Real-time model */
RT_MODEL_helicopter4_T helicopter4_M_;
RT_MODEL_helicopter4_T *const helicopter4_M = &helicopter4_M_;

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
  helicopter4_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter4_output(void)
{
  /* local block i/o variables */
  real_T rtb_FromWorkspace2[2];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Gain1_idx_4;
  real_T rtb_Gain1_idx_5;
  if (rtmIsMajorTimeStep(helicopter4_M)) {
    /* set solver stop time */
    if (!(helicopter4_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter4_M->solverInfo,
                            ((helicopter4_M->Timing.clockTickH0 + 1) *
        helicopter4_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter4_M->solverInfo,
                            ((helicopter4_M->Timing.clockTick0 + 1) *
        helicopter4_M->Timing.stepSize0 + helicopter4_M->Timing.clockTickH0 *
        helicopter4_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter4_M)) {
    helicopter4_M->Timing.t[0] = rtsiGetT(&helicopter4_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter4_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S6>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter4_DW.HILReadEncoderTimebase_Task,
        1, &helicopter4_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter4_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter4_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter4_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S6>/Travel: Count to rad' */
    helicopter4_B.TravelCounttorad = helicopter4_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S13>/Gain' */
    helicopter4_B.Gain = helicopter4_P.Gain_Gain *
      helicopter4_B.TravelCounttorad;

    /* Sum: '<Root>/Sum4' incorporates:
     *  Constant: '<Root>/ '
     */
    helicopter4_B.Sum4 = helicopter4_P._Value + helicopter4_B.Gain;

    /* Gain: '<S6>/Pitch: Count to rad' */
    helicopter4_B.PitchCounttorad = helicopter4_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S10>/Gain' */
    helicopter4_B.Gain_i = helicopter4_P.Gain_Gain_a *
      helicopter4_B.PitchCounttorad;

    /* Gain: '<S6>/Elevation: Count to rad' */
    helicopter4_B.ElevationCounttorad = helicopter4_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S8>/Gain' */
    helicopter4_B.Gain_e = helicopter4_P.Gain_Gain_l *
      helicopter4_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter4_B.Sum = helicopter4_B.Gain_e +
      helicopter4_P.elavation_offsetdeg_Value;

    /* Gain: '<S1>/Gain1' */
    helicopter4_B.Gain1[0] = helicopter4_P.Gain1_Gain * helicopter4_B.Sum4;
    helicopter4_B.Gain1[1] = helicopter4_P.Gain1_Gain * helicopter4_B.Gain_i;
    helicopter4_B.Gain1[2] = helicopter4_P.Gain1_Gain * helicopter4_B.Sum;
  }

  /* FromWorkspace: '<Root>/From Workspace2' */
  {
    real_T *pDataValues = (real_T *) helicopter4_DW.FromWorkspace2_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter4_DW.FromWorkspace2_PWORK.TimePtr;
    int_T currTimeIndex = helicopter4_DW.FromWorkspace2_IWORK.PrevIndex;
    real_T t = helicopter4_M->Timing.t[0];

    /* Get index */
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

    helicopter4_DW.FromWorkspace2_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace2[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace2[0])[elIdx] = pDataValues[currTimeIndex + 1];
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
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_FromWorkspace2[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter4_M)) {
    /* Gain: '<S4>/Gain1' incorporates:
     *  Constant: '<Root>/ 1'
     */
    helicopter4_B.Gain1_h = helicopter4_P.Gain1_Gain_l * helicopter4_P.u_Value;
  }

  /* TransferFcn: '<S6>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter4_P.TravelTransferFcn_C *
    helicopter4_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter4_P.TravelTransferFcn_D *
    helicopter4_B.TravelCounttorad;

  /* Gain: '<S14>/Gain' */
  helicopter4_B.Gain_d = helicopter4_P.Gain_Gain_lu * rtb_Backgain;

  /* TransferFcn: '<S6>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter4_P.PitchTransferFcn_C *
    helicopter4_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter4_P.PitchTransferFcn_D *
    helicopter4_B.PitchCounttorad;

  /* Gain: '<S11>/Gain' */
  helicopter4_B.Gain_b = helicopter4_P.Gain_Gain_ae * rtb_Backgain;

  /* TransferFcn: '<S6>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter4_P.ElevationTransferFcn_C *
    helicopter4_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter4_P.ElevationTransferFcn_D *
    helicopter4_B.ElevationCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter4_B.Gain_dg = helicopter4_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S3>/Gain1' */
  rtb_Gain1_idx_4 = helicopter4_P.Gain1_Gain_f * helicopter4_B.Sum;
  rtb_Gain1_idx_5 = helicopter4_P.Gain1_Gain_f * helicopter4_B.Gain_dg;

  /* Sum: '<S7>/Sum' incorporates:
   *  Constant: '<S7>/Vd_bias'
   *  Gain: '<S3>/Gain1'
   *  Gain: '<S7>/K_pd'
   *  Gain: '<S7>/K_pp'
   *  Sum: '<Root>/Sum1'
   *  Sum: '<S7>/Sum2'
   *  Sum: '<S7>/Sum3'
   */
  helicopter4_B.Sum_k = ((rtb_FromWorkspace2[0] - (helicopter4_P.Gain1_Gain_f *
    helicopter4_B.Gain_i + helicopter4_B.Gain1_h)) * helicopter4_P.K_pp -
    helicopter4_P.Gain1_Gain_f * helicopter4_B.Gain_b * helicopter4_P.K_pd) +
    helicopter4_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter4_M)) {
  }

  /* Integrator: '<S5>/Integrator'
   *
   * Regarding '<S5>/Integrator':
   *  Limited Integrator
   */
  if (helicopter4_X.Integrator_CSTATE >= helicopter4_P.Integrator_UpperSat ) {
    helicopter4_X.Integrator_CSTATE = helicopter4_P.Integrator_UpperSat;
  } else if (helicopter4_X.Integrator_CSTATE <=
             (helicopter4_P.Integrator_LowerSat) ) {
    helicopter4_X.Integrator_CSTATE = (helicopter4_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter4_X.Integrator_CSTATE;

  /* Sum: '<S5>/Sum' incorporates:
   *  Gain: '<Root>/Gain'
   */
  rtb_Gain1_idx_4 = helicopter4_P.Gain_Gain_d * rtb_FromWorkspace2[1] -
    rtb_Gain1_idx_4;

  /* Sum: '<S5>/Sum2' incorporates:
   *  Constant: '<S5>/Vs_bias'
   *  Gain: '<S5>/K_ed'
   *  Gain: '<S5>/K_ep'
   *  Sum: '<S5>/Sum1'
   */
  helicopter4_B.Sum2 = ((helicopter4_P.K_ep * rtb_Gain1_idx_4 + rtb_Backgain) -
                        helicopter4_P.K_ed * rtb_Gain1_idx_5) +
    helicopter4_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter4_M)) {
  }

  /* Gain: '<S2>/Back gain' incorporates:
   *  Sum: '<S2>/Subtract'
   */
  rtb_Backgain = (helicopter4_B.Sum2 - helicopter4_B.Sum_k) *
    helicopter4_P.Backgain_Gain;

  /* Gain: '<S5>/K_ei' */
  helicopter4_B.K_ei = helicopter4_P.K_ei * rtb_Gain1_idx_4;
  if (rtmIsMajorTimeStep(helicopter4_M)) {
  }

  /* Derivative: '<S6>/Derivative' */
  if ((helicopter4_DW.TimeStampA >= helicopter4_M->Timing.t[0]) &&
      (helicopter4_DW.TimeStampB >= helicopter4_M->Timing.t[0])) {
    rtb_Gain1_idx_4 = 0.0;
  } else {
    rtb_Gain1_idx_4 = helicopter4_DW.TimeStampA;
    lastU = &helicopter4_DW.LastUAtTimeA;
    if (helicopter4_DW.TimeStampA < helicopter4_DW.TimeStampB) {
      if (helicopter4_DW.TimeStampB < helicopter4_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helicopter4_DW.TimeStampB;
        lastU = &helicopter4_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter4_DW.TimeStampA >= helicopter4_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helicopter4_DW.TimeStampB;
        lastU = &helicopter4_DW.LastUAtTimeB;
      }
    }

    rtb_Gain1_idx_4 = (helicopter4_B.PitchCounttorad - *lastU) /
      (helicopter4_M->Timing.t[0] - rtb_Gain1_idx_4);
  }

  /* End of Derivative: '<S6>/Derivative' */

  /* Gain: '<S12>/Gain' */
  helicopter4_B.Gain_l = helicopter4_P.Gain_Gain_a1 * rtb_Gain1_idx_4;
  if (rtmIsMajorTimeStep(helicopter4_M)) {
  }

  /* Saturate: '<S6>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter4_P.BackmotorSaturation_UpperSat) {
    helicopter4_B.BackmotorSaturation =
      helicopter4_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter4_P.BackmotorSaturation_LowerSat) {
    helicopter4_B.BackmotorSaturation =
      helicopter4_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter4_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S6>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter4_M)) {
  }

  /* Gain: '<S2>/Front gain' incorporates:
   *  Sum: '<S2>/Add'
   */
  rtb_Gain1_idx_4 = (helicopter4_B.Sum_k + helicopter4_B.Sum2) *
    helicopter4_P.Frontgain_Gain;

  /* Saturate: '<S6>/Front motor: Saturation' */
  if (rtb_Gain1_idx_4 > helicopter4_P.FrontmotorSaturation_UpperSat) {
    helicopter4_B.FrontmotorSaturation =
      helicopter4_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 < helicopter4_P.FrontmotorSaturation_LowerSat) {
    helicopter4_B.FrontmotorSaturation =
      helicopter4_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter4_B.FrontmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S6>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter4_M)) {
    /* S-Function (hil_write_analog_block): '<S6>/HIL Write Analog' */

    /* S-Function Block: helicopter4/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter4_DW.HILWriteAnalog_Buffer[0] =
        helicopter4_B.FrontmotorSaturation;
      helicopter4_DW.HILWriteAnalog_Buffer[1] =
        helicopter4_B.BackmotorSaturation;
      result = hil_write_analog(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILWriteAnalog_channels, 2,
        &helicopter4_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter4_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S6>/Derivative' */
  if (helicopter4_DW.TimeStampA == (rtInf)) {
    helicopter4_DW.TimeStampA = helicopter4_M->Timing.t[0];
    lastU = &helicopter4_DW.LastUAtTimeA;
  } else if (helicopter4_DW.TimeStampB == (rtInf)) {
    helicopter4_DW.TimeStampB = helicopter4_M->Timing.t[0];
    lastU = &helicopter4_DW.LastUAtTimeB;
  } else if (helicopter4_DW.TimeStampA < helicopter4_DW.TimeStampB) {
    helicopter4_DW.TimeStampA = helicopter4_M->Timing.t[0];
    lastU = &helicopter4_DW.LastUAtTimeA;
  } else {
    helicopter4_DW.TimeStampB = helicopter4_M->Timing.t[0];
    lastU = &helicopter4_DW.LastUAtTimeB;
  }

  *lastU = helicopter4_B.PitchCounttorad;

  /* End of Update for Derivative: '<S6>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter4_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter4_M->solverInfo);
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
  if (!(++helicopter4_M->Timing.clockTick0)) {
    ++helicopter4_M->Timing.clockTickH0;
  }

  helicopter4_M->Timing.t[0] = rtsiGetSolverStopTime(&helicopter4_M->solverInfo);

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
    if (!(++helicopter4_M->Timing.clockTick1)) {
      ++helicopter4_M->Timing.clockTickH1;
    }

    helicopter4_M->Timing.t[1] = helicopter4_M->Timing.clockTick1 *
      helicopter4_M->Timing.stepSize1 + helicopter4_M->Timing.clockTickH1 *
      helicopter4_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter4_derivatives(void)
{
  XDot_helicopter4_T *_rtXdot;
  _rtXdot = ((XDot_helicopter4_T *) helicopter4_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S6>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter4_P.TravelTransferFcn_A *
    helicopter4_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter4_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S6>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter4_P.PitchTransferFcn_A *
    helicopter4_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter4_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S6>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter4_P.ElevationTransferFcn_A *
    helicopter4_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter4_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S5>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter4_X.Integrator_CSTATE <=
            (helicopter4_P.Integrator_LowerSat) );
    usat = ( helicopter4_X.Integrator_CSTATE >=
            helicopter4_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter4_B.K_ei > 0)) ||
        (usat && (helicopter4_B.K_ei < 0)) ) {
      ((XDot_helicopter4_T *) helicopter4_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter4_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter4_T *) helicopter4_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter4_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter4/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter4_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter4_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter4_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter4_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter4_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter4_M, _rt_error_message);
      return;
    }

    if ((helicopter4_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter4_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter4_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter4_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter4_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter4_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_analog_input_chan, 8U,
        &helicopter4_DW.HILInitialize_AIMinimums[0],
        &helicopter4_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter4_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter4_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter4_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter4_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter4_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter4_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_analog_output_cha, 8U,
        &helicopter4_DW.HILInitialize_AOMinimums[0],
        &helicopter4_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter4_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter4_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter4_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_analog_output_cha, 8U,
        &helicopter4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if (helicopter4_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter4_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter4_DW.HILInitialize_Card,
         helicopter4_P.HILInitialize_analog_output_cha, 8U,
         &helicopter4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter4_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter4_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter4_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter4_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_encoder_channels, 8U,
        (t_encoder_quadrature_mode *)
        &helicopter4_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter4_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter4_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter4_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicopter4_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_encoder_channels, 8U,
        &helicopter4_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter4_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter4_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter4_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter4_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter4_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helicopter4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter4_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter4_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter4_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helicopter4_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter4_DW.HILInitialize_Card,
          &helicopter4_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter4_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter4_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter4_DW.HILInitialize_Card,
          &helicopter4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter4_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter4_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &helicopter4_DW.HILInitialize_POAlignValues
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter4_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter4_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter4_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter4_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter4_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter4_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter4_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter4_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter4_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_pwm_channels, 8U,
        &helicopter4_DW.HILInitialize_POSortedFreqs[0],
        &helicopter4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter4_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter4_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter4_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter4_DW.HILInitialize_Card,
        helicopter4_P.HILInitialize_pwm_channels, 8U,
        &helicopter4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }

    if (helicopter4_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter4_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter4_DW.HILInitialize_Card,
         helicopter4_P.HILInitialize_pwm_channels, 8U,
         &helicopter4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S6>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter4_DW.HILInitialize_Card,
      helicopter4_P.HILReadEncoderTimebase_samples_,
      helicopter4_P.HILReadEncoderTimebase_channels, 3,
      &helicopter4_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter4_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace2' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359867411190242,
      0.52359865554758189, 0.52359863151719044, 0.52359859964940625,
      0.52359855608761319, 0.52359849415401238, 0.523598401303277,
      0.5235982505930511, 0.52359797202561931, 0.52359731103603413,
      0.52359399929140371, 0.37373278194791881, 0.10884456730038469,
      -0.10901743515278248, -0.2741585625437295, -0.39337971917392456,
      -0.47327537098995542, -0.51994722019507711, -0.52354978765231674,
      -0.523542767065935, -0.5186834872817716, -0.48724613927396954,
      -0.44678840376950552, -0.40067266669214524, -0.35159284154915316,
      -0.30178334627213843, -0.25305064324654308, -0.20681491122926998,
      -0.16417308677197467, -0.12592972817989473, -0.092639633165663912,
      -0.064634795410435833, -0.042051893964153819, -0.024828919455665236,
      -0.012683447163317594, -0.0050915690677124864, -0.0012398265018539612,
      1.2200584620951809E-6, -1.1972103431244371E-6, 2.9754361053828935E-6,
      2.9754361053828935E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.28557759301561392, 0.30305180452654246, 0.31969613804665942,
      0.33492015191161323, 0.34796549534863969, 0.35790622068453448,
      0.36358792731773448, 0.36358805659499566, 0.35616764841922405,
      0.33920922806101789, 0.3101347272678962, 0.26583073130615309,
      0.20254119862726033, 0.11574263916140465, 1.2613648081217309E-5,
      1.5697396982913826E-5, 1.4460416593932361E-5, 1.5310600547693464E-5,
      1.4787166619899529E-5, 1.4643583279962577E-5, 1.7012691720968447E-5,
      1.4768551726623247E-5, 1.426590543664584E-5, 1.5614096264903502E-5,
      1.5891102263351359E-5, 1.6098563914480772E-5, 1.5390740083914525E-5,
      1.3674933093567232E-5, 1.5122495924767908E-5, 1.164780967720294E-5,
      1.1232280046905501E-5, 1.0123615178452864E-5, 7.9133678004933361E-6,
      7.8082546513002624E-6, 6.1501733234757783E-6, 4.0454039505571296E-6,
      3.9033798065508141E-6, 1.0150832640510402E-6, -4.1535009460853133E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter4_DW.FromWorkspace2_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter4_DW.FromWorkspace2_PWORK.DataPtr = (void *) pDataValues0;
    helicopter4_DW.FromWorkspace2_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S6>/Travel: Transfer Fcn' */
  helicopter4_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S6>/Pitch: Transfer Fcn' */
  helicopter4_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S6>/Elevation: Transfer Fcn' */
  helicopter4_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S5>/Integrator' */
  helicopter4_X.Integrator_CSTATE = helicopter4_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S6>/Derivative' */
  helicopter4_DW.TimeStampA = (rtInf);
  helicopter4_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter4_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter4/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter4_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter4_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter4_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter4_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter4_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter4_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter4_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter4_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter4_DW.HILInitialize_Card
                         , helicopter4_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter4_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter4_DW.HILInitialize_AOVoltages[0]
                         , &helicopter4_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter4_DW.HILInitialize_Card,
            helicopter4_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helicopter4_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter4_DW.HILInitialize_Card,
            helicopter4_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter4_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter4_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter4_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter4_DW.HILInitialize_Card);
    hil_close(helicopter4_DW.HILInitialize_Card);
    helicopter4_DW.HILInitialize_Card = NULL;
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
  helicopter4_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter4_update();
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
  helicopter4_initialize();
}

void MdlTerminate(void)
{
  helicopter4_terminate();
}

/* Registration function */
RT_MODEL_helicopter4_T *helicopter4(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter4_P.Integrator_UpperSat = rtInf;
  helicopter4_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter4_M, 0,
                sizeof(RT_MODEL_helicopter4_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter4_M->solverInfo,
                          &helicopter4_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter4_M->solverInfo, &rtmGetTPtr(helicopter4_M));
    rtsiSetStepSizePtr(&helicopter4_M->solverInfo,
                       &helicopter4_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter4_M->solverInfo, &helicopter4_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter4_M->solverInfo, (real_T **)
                         &helicopter4_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter4_M->solverInfo,
      &helicopter4_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter4_M->solverInfo, (&rtmGetErrorStatus
      (helicopter4_M)));
    rtsiSetRTModelPtr(&helicopter4_M->solverInfo, helicopter4_M);
  }

  rtsiSetSimTimeStep(&helicopter4_M->solverInfo, MAJOR_TIME_STEP);
  helicopter4_M->ModelData.intgData.f[0] = helicopter4_M->ModelData.odeF[0];
  helicopter4_M->ModelData.contStates = ((real_T *) &helicopter4_X);
  rtsiSetSolverData(&helicopter4_M->solverInfo, (void *)
                    &helicopter4_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter4_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter4_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter4_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter4_M->Timing.sampleTimes = (&helicopter4_M->
      Timing.sampleTimesArray[0]);
    helicopter4_M->Timing.offsetTimes = (&helicopter4_M->
      Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter4_M->Timing.sampleTimes[0] = (0.0);
    helicopter4_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter4_M->Timing.offsetTimes[0] = (0.0);
    helicopter4_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter4_M, &helicopter4_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter4_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter4_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter4_M, 20.0);
  helicopter4_M->Timing.stepSize0 = 0.002;
  helicopter4_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter4_M->Sizes.checksums[0] = (3419306624U);
  helicopter4_M->Sizes.checksums[1] = (2791415762U);
  helicopter4_M->Sizes.checksums[2] = (3013053323U);
  helicopter4_M->Sizes.checksums[3] = (2715059400U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter4_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter4_M->extModeInfo,
      &helicopter4_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter4_M->extModeInfo,
                        helicopter4_M->Sizes.checksums);
    rteiSetTPtr(helicopter4_M->extModeInfo, rtmGetTPtr(helicopter4_M));
  }

  helicopter4_M->solverInfoPtr = (&helicopter4_M->solverInfo);
  helicopter4_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter4_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter4_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter4_M->ModelData.blockIO = ((void *) &helicopter4_B);

  {
    helicopter4_B.TravelCounttorad = 0.0;
    helicopter4_B.Gain = 0.0;
    helicopter4_B.Sum4 = 0.0;
    helicopter4_B.PitchCounttorad = 0.0;
    helicopter4_B.Gain_i = 0.0;
    helicopter4_B.ElevationCounttorad = 0.0;
    helicopter4_B.Gain_e = 0.0;
    helicopter4_B.Sum = 0.0;
    helicopter4_B.Gain1[0] = 0.0;
    helicopter4_B.Gain1[1] = 0.0;
    helicopter4_B.Gain1[2] = 0.0;
    helicopter4_B.Gain1_h = 0.0;
    helicopter4_B.Gain_d = 0.0;
    helicopter4_B.Gain_b = 0.0;
    helicopter4_B.Gain_dg = 0.0;
    helicopter4_B.Sum_k = 0.0;
    helicopter4_B.Sum2 = 0.0;
    helicopter4_B.K_ei = 0.0;
    helicopter4_B.Gain_l = 0.0;
    helicopter4_B.BackmotorSaturation = 0.0;
    helicopter4_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter4_M->ModelData.defaultParam = ((real_T *)&helicopter4_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter4_X;
    helicopter4_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter4_X, 0,
                  sizeof(X_helicopter4_T));
  }

  /* states (dwork) */
  helicopter4_M->ModelData.dwork = ((void *) &helicopter4_DW);
  (void) memset((void *)&helicopter4_DW, 0,
                sizeof(DW_helicopter4_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter4_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter4_DW.TimeStampA = 0.0;
  helicopter4_DW.LastUAtTimeA = 0.0;
  helicopter4_DW.TimeStampB = 0.0;
  helicopter4_DW.LastUAtTimeB = 0.0;
  helicopter4_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter4_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter4_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter4_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter4_M->Sizes.numY = (0);     /* Number of model outputs */
  helicopter4_M->Sizes.numU = (0);     /* Number of model inputs */
  helicopter4_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter4_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter4_M->Sizes.numBlocks = (60);/* Number of blocks */
  helicopter4_M->Sizes.numBlockIO = (19);/* Number of block outputs */
  helicopter4_M->Sizes.numBlockPrms = (145);/* Sum of parameter "widths" */
  return helicopter4_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
