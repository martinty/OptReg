/*
 * helicopter44.c
 *
 * Code generation for model "helicopter44".
 *
 * Model version              : 1.185
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Thu Apr 27 15:14:37 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter44.h"
#include "helicopter44_private.h"
#include "helicopter44_dt.h"

/* Block signals (auto storage) */
B_helicopter44_T helicopter44_B;

/* Continuous states */
X_helicopter44_T helicopter44_X;

/* Block states (auto storage) */
DW_helicopter44_T helicopter44_DW;

/* Real-time model */
RT_MODEL_helicopter44_T helicopter44_M_;
RT_MODEL_helicopter44_T *const helicopter44_M = &helicopter44_M_;

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
  helicopter44_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter44_output(void)
{
  /* local block i/o variables */
  real_T rtb_FromWorkspace2[2];
  real_T rtb_Sum1_d[6];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Gain1[6];
  real_T rtb_Sum3_m[2];
  real_T rtb_Derivative;
  int32_T i;
  int32_T i_0;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* set solver stop time */
    if (!(helicopter44_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter44_M->solverInfo,
                            ((helicopter44_M->Timing.clockTickH0 + 1) *
        helicopter44_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter44_M->solverInfo,
                            ((helicopter44_M->Timing.clockTick0 + 1) *
        helicopter44_M->Timing.stepSize0 + helicopter44_M->Timing.clockTickH0 *
        helicopter44_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter44_M)) {
    helicopter44_M->Timing.t[0] = rtsiGetT(&helicopter44_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter44/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter44_DW.HILReadEncoderTimebase_Task,
        1, &helicopter44_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter44_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter44_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter44_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S5>/Travel: Count to rad' */
    helicopter44_B.TravelCounttorad = helicopter44_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helicopter44_B.Gain = helicopter44_P.Gain_Gain *
      helicopter44_B.TravelCounttorad;

    /* Sum: '<Root>/Sum4' incorporates:
     *  Constant: '<Root>/ '
     */
    helicopter44_B.Sum4 = helicopter44_P._Value + helicopter44_B.Gain;

    /* Gain: '<S5>/Pitch: Count to rad' */
    helicopter44_B.PitchCounttorad = helicopter44_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helicopter44_B.Gain_i = helicopter44_P.Gain_Gain_a *
      helicopter44_B.PitchCounttorad;

    /* Gain: '<S5>/Elevation: Count to rad' */
    helicopter44_B.ElevationCounttorad = helicopter44_P.ElevationCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    helicopter44_B.Gain_e = helicopter44_P.Gain_Gain_l *
      helicopter44_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter44_B.Sum = helicopter44_B.Gain_e +
      helicopter44_P.elavation_offsetdeg_Value;

    /* Gain: '<S1>/Gain1' */
    helicopter44_B.Gain1[0] = helicopter44_P.Gain1_Gain * helicopter44_B.Sum4;
    helicopter44_B.Gain1[1] = helicopter44_P.Gain1_Gain * helicopter44_B.Gain_i;
    helicopter44_B.Gain1[2] = helicopter44_P.Gain1_Gain * helicopter44_B.Sum;
  }

  /* FromWorkspace: '<Root>/From Workspace3' */
  {
    real_T *pDataValues = (real_T *)
      helicopter44_DW.FromWorkspace3_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter44_DW.FromWorkspace3_PWORK.TimePtr;
    int_T currTimeIndex = helicopter44_DW.FromWorkspace3_IWORK.PrevIndex;
    real_T t = helicopter44_M->Timing.t[0];

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

    helicopter44_DW.FromWorkspace3_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum1_d[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum1_d[0])[elIdx] = pDataValues[currTimeIndex + 1];
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
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum1_d[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  /* TransferFcn: '<S5>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter44_P.TravelTransferFcn_C *
    helicopter44_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter44_P.TravelTransferFcn_D *
    helicopter44_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  helicopter44_B.Gain_d = helicopter44_P.Gain_Gain_lu * rtb_Backgain;

  /* TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter44_P.PitchTransferFcn_C *
    helicopter44_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter44_P.PitchTransferFcn_D *
    helicopter44_B.PitchCounttorad;

  /* Gain: '<S10>/Gain' */
  helicopter44_B.Gain_b = helicopter44_P.Gain_Gain_ae * rtb_Backgain;

  /* TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter44_P.ElevationTransferFcn_C *
    helicopter44_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter44_P.ElevationTransferFcn_D *
    helicopter44_B.ElevationCounttorad;

  /* Gain: '<S8>/Gain' */
  helicopter44_B.Gain_dg = helicopter44_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S3>/Gain1' */
  rtb_Gain1[0] = helicopter44_P.Gain1_Gain_f * helicopter44_B.Sum4;
  rtb_Gain1[1] = helicopter44_P.Gain1_Gain_f * helicopter44_B.Gain_d;
  rtb_Gain1[2] = helicopter44_P.Gain1_Gain_f * helicopter44_B.Gain_i;
  rtb_Gain1[3] = helicopter44_P.Gain1_Gain_f * helicopter44_B.Gain_b;
  rtb_Gain1[4] = helicopter44_P.Gain1_Gain_f * helicopter44_B.Sum;
  rtb_Gain1[5] = helicopter44_P.Gain1_Gain_f * helicopter44_B.Gain_dg;

  /* Sum: '<Root>/Sum1' */
  for (i = 0; i < 6; i++) {
    rtb_Sum1_d[i] -= rtb_Gain1[i];
  }

  /* End of Sum: '<Root>/Sum1' */

  /* FromWorkspace: '<Root>/From Workspace2' */
  {
    real_T *pDataValues = (real_T *)
      helicopter44_DW.FromWorkspace2_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter44_DW.FromWorkspace2_PWORK.TimePtr;
    int_T currTimeIndex = helicopter44_DW.FromWorkspace2_IWORK.PrevIndex;
    real_T t = helicopter44_M->Timing.t[0];

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

    helicopter44_DW.FromWorkspace2_IWORK.PrevIndex = currTimeIndex;

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

  /* Sum: '<Root>/Sum3' incorporates:
   *  Gain: '<Root>/Gain'
   */
  for (i = 0; i < 2; i++) {
    rtb_Derivative = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Derivative += helicopter44_P.K_LQ[(i_0 << 1) + i] * rtb_Sum1_d[i_0];
    }

    rtb_Sum3_m[i] = rtb_Derivative + rtb_FromWorkspace2[i];
  }

  /* End of Sum: '<Root>/Sum3' */

  /* Sum: '<S6>/Sum' incorporates:
   *  Constant: '<S6>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  helicopter44_B.Sum_k = ((rtb_Sum3_m[0] - rtb_Gain1[2]) * helicopter44_P.K_pp -
    helicopter44_P.K_pd * rtb_Gain1[3]) + helicopter44_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Integrator: '<S4>/Integrator'
   *
   * Regarding '<S4>/Integrator':
   *  Limited Integrator
   */
  if (helicopter44_X.Integrator_CSTATE >= helicopter44_P.Integrator_UpperSat ) {
    helicopter44_X.Integrator_CSTATE = helicopter44_P.Integrator_UpperSat;
  } else if (helicopter44_X.Integrator_CSTATE <=
             (helicopter44_P.Integrator_LowerSat) ) {
    helicopter44_X.Integrator_CSTATE = (helicopter44_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter44_X.Integrator_CSTATE;

  /* Sum: '<S4>/Sum' incorporates:
   *  Gain: '<Root>/        '
   */
  rtb_Derivative = helicopter44_P._Gain * rtb_Sum3_m[1] - rtb_Gain1[4];

  /* Sum: '<S4>/Sum2' incorporates:
   *  Constant: '<S4>/Vs_bias'
   *  Gain: '<S4>/K_ed'
   *  Gain: '<S4>/K_ep'
   *  Sum: '<S4>/Sum1'
   */
  helicopter44_B.Sum2 = ((helicopter44_P.K_ep * rtb_Derivative + rtb_Backgain) -
    helicopter44_P.K_ed * rtb_Gain1[5]) + helicopter44_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Gain: '<S2>/Back gain' incorporates:
   *  Sum: '<S2>/Subtract'
   */
  rtb_Backgain = (helicopter44_B.Sum2 - helicopter44_B.Sum_k) *
    helicopter44_P.Backgain_Gain;

  /* Gain: '<S4>/K_ei' */
  helicopter44_B.K_ei = helicopter44_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Derivative: '<S5>/Derivative' */
  if ((helicopter44_DW.TimeStampA >= helicopter44_M->Timing.t[0]) &&
      (helicopter44_DW.TimeStampB >= helicopter44_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter44_DW.TimeStampA;
    lastU = &helicopter44_DW.LastUAtTimeA;
    if (helicopter44_DW.TimeStampA < helicopter44_DW.TimeStampB) {
      if (helicopter44_DW.TimeStampB < helicopter44_M->Timing.t[0]) {
        rtb_Derivative = helicopter44_DW.TimeStampB;
        lastU = &helicopter44_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter44_DW.TimeStampA >= helicopter44_M->Timing.t[0]) {
        rtb_Derivative = helicopter44_DW.TimeStampB;
        lastU = &helicopter44_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter44_B.PitchCounttorad - *lastU) /
      (helicopter44_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S5>/Derivative' */

  /* Gain: '<S11>/Gain' */
  helicopter44_B.Gain_l = helicopter44_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Saturate: '<S5>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter44_P.BackmotorSaturation_UpperSat) {
    helicopter44_B.BackmotorSaturation =
      helicopter44_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter44_P.BackmotorSaturation_LowerSat) {
    helicopter44_B.BackmotorSaturation =
      helicopter44_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter44_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S5>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Gain: '<S2>/Front gain' incorporates:
   *  Sum: '<S2>/Add'
   */
  rtb_Derivative = (helicopter44_B.Sum_k + helicopter44_B.Sum2) *
    helicopter44_P.Frontgain_Gain;

  /* Saturate: '<S5>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter44_P.FrontmotorSaturation_UpperSat) {
    helicopter44_B.FrontmotorSaturation =
      helicopter44_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter44_P.FrontmotorSaturation_LowerSat) {
    helicopter44_B.FrontmotorSaturation =
      helicopter44_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter44_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S5>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* S-Function (hil_write_analog_block): '<S5>/HIL Write Analog' */

    /* S-Function Block: helicopter44/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter44_DW.HILWriteAnalog_Buffer[0] =
        helicopter44_B.FrontmotorSaturation;
      helicopter44_DW.HILWriteAnalog_Buffer[1] =
        helicopter44_B.BackmotorSaturation;
      result = hil_write_analog(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILWriteAnalog_channels, 2,
        &helicopter44_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter44_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S5>/Derivative' */
  if (helicopter44_DW.TimeStampA == (rtInf)) {
    helicopter44_DW.TimeStampA = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeA;
  } else if (helicopter44_DW.TimeStampB == (rtInf)) {
    helicopter44_DW.TimeStampB = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeB;
  } else if (helicopter44_DW.TimeStampA < helicopter44_DW.TimeStampB) {
    helicopter44_DW.TimeStampA = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeA;
  } else {
    helicopter44_DW.TimeStampB = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeB;
  }

  *lastU = helicopter44_B.PitchCounttorad;

  /* End of Update for Derivative: '<S5>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter44_M->solverInfo);
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
  if (!(++helicopter44_M->Timing.clockTick0)) {
    ++helicopter44_M->Timing.clockTickH0;
  }

  helicopter44_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter44_M->solverInfo);

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
    if (!(++helicopter44_M->Timing.clockTick1)) {
      ++helicopter44_M->Timing.clockTickH1;
    }

    helicopter44_M->Timing.t[1] = helicopter44_M->Timing.clockTick1 *
      helicopter44_M->Timing.stepSize1 + helicopter44_M->Timing.clockTickH1 *
      helicopter44_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter44_derivatives(void)
{
  XDot_helicopter44_T *_rtXdot;
  _rtXdot = ((XDot_helicopter44_T *) helicopter44_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter44_P.TravelTransferFcn_A *
    helicopter44_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter44_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter44_P.PitchTransferFcn_A *
    helicopter44_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter44_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter44_P.ElevationTransferFcn_A *
    helicopter44_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter44_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S4>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter44_X.Integrator_CSTATE <=
            (helicopter44_P.Integrator_LowerSat) );
    usat = ( helicopter44_X.Integrator_CSTATE >=
            helicopter44_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter44_B.K_ei > 0)) ||
        (usat && (helicopter44_B.K_ei < 0)) ) {
      ((XDot_helicopter44_T *) helicopter44_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter44_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter44_T *) helicopter44_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter44_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter44/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter44_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter44_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter44_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      return;
    }

    if ((helicopter44_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter44_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter44_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter44_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter44_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_analog_input_chan, 8U,
        &helicopter44_DW.HILInitialize_AIMinimums[0],
        &helicopter44_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter44_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter44_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter44_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter44_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_analog_output_cha, 8U,
        &helicopter44_DW.HILInitialize_AOMinimums[0],
        &helicopter44_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter44_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter44_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_analog_output_cha, 8U,
        &helicopter44_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if (helicopter44_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter44_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter44_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter44_DW.HILInitialize_Card,
         helicopter44_P.HILInitialize_analog_output_cha, 8U,
         &helicopter44_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter44_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter44_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter44_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter44_DW.HILInitialize_Card,
         helicopter44_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter44_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter44_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter44_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter44_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_encoder_channels, 8U,
        &helicopter44_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter44_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter44_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter44_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter44_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter44_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helicopter44_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter44_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter44_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter44_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter44_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter44_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helicopter44_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter44_DW.HILInitialize_Card,
          &helicopter44_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter44_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter44_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter44_DW.HILInitialize_Card,
          &helicopter44_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter44_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter44_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter44_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter44_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter44_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter44_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter44_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter44_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter44_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter44_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter44_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter44_DW.HILInitialize_POSortedFreqs
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter44_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U,
        &helicopter44_DW.HILInitialize_POSortedFreqs[0],
        &helicopter44_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter44_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U,
        &helicopter44_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if (helicopter44_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter44_DW.HILInitialize_Card,
         helicopter44_P.HILInitialize_pwm_channels, 8U,
         &helicopter44_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter44/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter44_DW.HILInitialize_Card,
      helicopter44_P.HILReadEncoderTimebase_samples_,
      helicopter44_P.HILReadEncoderTimebase_channels, 3,
      &helicopter44_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace3' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.14159265357763, 3.1415926537670646,
      3.1415926538719172, 3.1378421434505541, 3.126215557612666,
      3.10330931007741, 3.066627423683447, 3.0144539279117768,
      2.9456562386258036, 2.85950765598972, 2.7555512956394872,
      2.6335045746505492, 2.4931946209899212, 2.3345171181132405,
      2.1584826519823448, 1.9681323624724858, 1.7680657423335424,
      1.5636077548816543, 1.3600737077137603, 1.1622752739025943,
      0.97423054318298652, 0.79893023268390639, 0.6384481991056391,
      0.49412651727848428, 0.36658871869139026, 0.2558250738787905,
      0.16130583596972617, 0.082099731413308757, 0.01698578810762931,
      -0.035446404322311519, -0.076710147751138755, -0.10834886177550877,
      -0.13187647153807305, -0.14872986466854082, -0.16023225278329037,
      -0.16756627024450593, -0.17175591068223625, -0.17365681760988025,
      -0.17395482106770829, -0.173173091779404, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 2.3913094913698171E-10, 8.002818257358723E-10,
      -0.01500204533746718, -0.046506341016376178, -0.0916249980861316,
      -0.14672751552647986, -0.20869395885386088, -0.275190575929402,
      -0.34459411173452947, -0.41582478847932852, -0.48818611621229135,
      -0.561238204548677, -0.63470828035649351, -0.704134959634044,
      -0.761398356705758, -0.80026250597829518, -0.81783584531711484,
      -0.81413894148916122, -0.79119657837785173, -0.75218066356526558,
      -0.7012028334383783, -0.64192895245670245, -0.57728739257606032,
      -0.51015146701021086, -0.44305477840211843, -0.37807700943063705,
      -0.31682445819252092, -0.26045577873506554, -0.20972877476271518,
      -0.16505497339955746, -0.12655485660126289, -0.094110438971837226,
      -0.067413572586396284, -0.046009552459765062, -0.029336069844861097,
      -0.016758561750921878, -0.0076036277105764232, -0.0011920138313117766,
      0.0031269171532173259, 0.00593819308084707, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -1.355030956579588E-11, 0.10602873111902804,
      0.22266033278758904, 0.31888139604824639, 0.3894435071209057,
      0.43795493670440028, 0.46997248881025339, 0.49051702264376218,
      0.50343076198511694, 0.511420954163112, 0.516303839459644,
      0.51925698244280682, 0.49068187165714233, 0.40471550959931113,
      0.27468004714709165, 0.1242041102733857, -0.026125208185550964,
      -0.16214803363454316, -0.27574894980004772, -0.36029311048436441,
      -0.4189255882504333, -0.45686301224152864, -0.47449158901540206,
      -0.47421439292458878, -0.45923831765367967, -0.43290986275339677,
      -0.39839237992174592, -0.358519143152337, -0.31573738104568055,
      -0.27210412885177226, -0.22930475379324142, -0.18868325593618585,
      -0.15127543956840997, -0.11784180713218065, -0.088893023565727189,
      -0.064703577315386743, -0.045314838154504988, -0.030524554698676268,
      -0.019869024564967658, -0.012622491660209286, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.424114926031222, 0.46652640359828929,
      0.38488425954066513, 0.28224842264941852, 0.19404574152066315,
      0.1280701175089419, 0.082178216138984669, 0.051654703228877982,
      0.0319610278746409, 0.019531044257214074, 0.011813227778291469,
      -0.11430116357533941, -0.34386427572733391, -0.52014275267570742,
      -0.60190243401010612, -0.60131849349105948, -0.54409040492116456,
      -0.45440470823866058, -0.33817604270492646, -0.23453047649278205,
      -0.15174934675030341, -0.070514512051594369, 0.0011089413253883539,
      0.059904252242653687, 0.1053138704194698, 0.13806992438019827,
      0.15949295831151586, 0.17112704799898371, 0.17453301032485022,
      0.17119750024115057, 0.16248599153176654, 0.14963126547231329,
      0.1337345297449167, 0.11579513426581486, 0.096757785001360688,
      0.077554956643527781, 0.059161133823314491, 0.042622120534834711,
      0.028986131619033368, 0.018994933162772461, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.1032902666612636E-10, 0.0044621502992741339,
      0.012543948046526794, 0.023530829209072839, 0.036808123159838849,
      0.051835397464152104, 0.068123040425041526, 0.08520997041635868,
      0.10264194026589683, 0.11994987189456932, 0.13662759060455873,
      0.15210813213633537, 0.16573819474794518, 0.17674981562118044,
      0.18422855383080416, 0.18707721532417754, 0.18633634965251014,
      0.18285246085319351, 0.17732775146868912, 0.17032709499706641,
      0.16230608302203772, 0.15362935616326859, 0.14458620952937334,
      0.13540380388250103, 0.12625823923418938, 0.1172837341077102,
      0.10858038520815894, 0.10022058838740507, 0.092254401261237656,
      0.084714056501418056, 0.077617507073612146, 0.070971613996005115,
      0.064774579032167129, 0.059017995036883361, 0.053688576270584022,
      0.048769452123524661, 0.044241288218610189, 0.040083203590885158,
      0.036273385852418953, 0.032789716002644784, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.017848599902222248, 0.032327188617228761,
      0.043947518316153039, 0.053109158465309982, 0.060109052967439208,
      0.065150461205291332, 0.068347465070098784, 0.0697273295768347,
      0.069230648954839408, 0.06670891936663248, 0.061918936794409976,
      0.054515297678667543, 0.044039561275136377, 0.029905900526938606,
      0.011383955893424646, -0.0029754608446538231, -0.013924242495596808,
      -0.022089044292642742, -0.027994490369295379, -0.0320779394304612,
      -0.034702608987612842, -0.036169841336661207, -0.036728001279519813,
      -0.036581398146334917, -0.035897604479608533, -0.034813218884075146,
      -0.033439119670701366, -0.031864726548167091, -0.030161372246036912,
      -0.028386195794529633, -0.026583571526281241, -0.024788139593078719,
      -0.023026335861399069, -0.021317675048440534, -0.019676496588322479,
      -0.018112655620169931, -0.016632338511435252, -0.01523927095439776,
      -0.01393467939952109, -0.01271809616538388, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter44_DW.FromWorkspace3_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter44_DW.FromWorkspace3_PWORK.DataPtr = (void *) pDataValues0;
    helicopter44_DW.FromWorkspace3_IWORK.PrevIndex = 0;
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

    helicopter44_DW.FromWorkspace2_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter44_DW.FromWorkspace2_PWORK.DataPtr = (void *) pDataValues0;
    helicopter44_DW.FromWorkspace2_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  helicopter44_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  helicopter44_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  helicopter44_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S4>/Integrator' */
  helicopter44_X.Integrator_CSTATE = helicopter44_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S5>/Derivative' */
  helicopter44_DW.TimeStampA = (rtInf);
  helicopter44_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter44_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter44/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter44_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter44_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter44_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter44_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter44_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter44_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter44_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter44_DW.HILInitialize_Card
                         , helicopter44_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter44_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter44_DW.HILInitialize_AOVoltages[0]
                         , &helicopter44_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter44_DW.HILInitialize_Card,
            helicopter44_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helicopter44_DW.HILInitialize_AOVoltages
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter44_DW.HILInitialize_Card,
            helicopter44_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter44_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter44_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter44_DW.HILInitialize_Card);
    hil_close(helicopter44_DW.HILInitialize_Card);
    helicopter44_DW.HILInitialize_Card = NULL;
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
  helicopter44_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter44_update();
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
  helicopter44_initialize();
}

void MdlTerminate(void)
{
  helicopter44_terminate();
}

/* Registration function */
RT_MODEL_helicopter44_T *helicopter44(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter44_P.Integrator_UpperSat = rtInf;
  helicopter44_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter44_M, 0,
                sizeof(RT_MODEL_helicopter44_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter44_M->solverInfo,
                          &helicopter44_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter44_M->solverInfo, &rtmGetTPtr(helicopter44_M));
    rtsiSetStepSizePtr(&helicopter44_M->solverInfo,
                       &helicopter44_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter44_M->solverInfo, &helicopter44_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter44_M->solverInfo, (real_T **)
                         &helicopter44_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter44_M->solverInfo,
      &helicopter44_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter44_M->solverInfo, (&rtmGetErrorStatus
      (helicopter44_M)));
    rtsiSetRTModelPtr(&helicopter44_M->solverInfo, helicopter44_M);
  }

  rtsiSetSimTimeStep(&helicopter44_M->solverInfo, MAJOR_TIME_STEP);
  helicopter44_M->ModelData.intgData.f[0] = helicopter44_M->ModelData.odeF[0];
  helicopter44_M->ModelData.contStates = ((real_T *) &helicopter44_X);
  rtsiSetSolverData(&helicopter44_M->solverInfo, (void *)
                    &helicopter44_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter44_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter44_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter44_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter44_M->Timing.sampleTimes =
      (&helicopter44_M->Timing.sampleTimesArray[0]);
    helicopter44_M->Timing.offsetTimes =
      (&helicopter44_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter44_M->Timing.sampleTimes[0] = (0.0);
    helicopter44_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter44_M->Timing.offsetTimes[0] = (0.0);
    helicopter44_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter44_M, &helicopter44_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter44_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter44_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter44_M, 20.0);
  helicopter44_M->Timing.stepSize0 = 0.002;
  helicopter44_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter44_M->Sizes.checksums[0] = (1449653811U);
  helicopter44_M->Sizes.checksums[1] = (2269923027U);
  helicopter44_M->Sizes.checksums[2] = (677359681U);
  helicopter44_M->Sizes.checksums[3] = (1393314079U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter44_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter44_M->extModeInfo,
      &helicopter44_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter44_M->extModeInfo,
                        helicopter44_M->Sizes.checksums);
    rteiSetTPtr(helicopter44_M->extModeInfo, rtmGetTPtr(helicopter44_M));
  }

  helicopter44_M->solverInfoPtr = (&helicopter44_M->solverInfo);
  helicopter44_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter44_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter44_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter44_M->ModelData.blockIO = ((void *) &helicopter44_B);

  {
    helicopter44_B.TravelCounttorad = 0.0;
    helicopter44_B.Gain = 0.0;
    helicopter44_B.Sum4 = 0.0;
    helicopter44_B.PitchCounttorad = 0.0;
    helicopter44_B.Gain_i = 0.0;
    helicopter44_B.ElevationCounttorad = 0.0;
    helicopter44_B.Gain_e = 0.0;
    helicopter44_B.Sum = 0.0;
    helicopter44_B.Gain1[0] = 0.0;
    helicopter44_B.Gain1[1] = 0.0;
    helicopter44_B.Gain1[2] = 0.0;
    helicopter44_B.Gain_d = 0.0;
    helicopter44_B.Gain_b = 0.0;
    helicopter44_B.Gain_dg = 0.0;
    helicopter44_B.Sum_k = 0.0;
    helicopter44_B.Sum2 = 0.0;
    helicopter44_B.K_ei = 0.0;
    helicopter44_B.Gain_l = 0.0;
    helicopter44_B.BackmotorSaturation = 0.0;
    helicopter44_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter44_M->ModelData.defaultParam = ((real_T *)&helicopter44_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter44_X;
    helicopter44_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter44_X, 0,
                  sizeof(X_helicopter44_T));
  }

  /* states (dwork) */
  helicopter44_M->ModelData.dwork = ((void *) &helicopter44_DW);
  (void) memset((void *)&helicopter44_DW, 0,
                sizeof(DW_helicopter44_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter44_DW.TimeStampA = 0.0;
  helicopter44_DW.LastUAtTimeA = 0.0;
  helicopter44_DW.TimeStampB = 0.0;
  helicopter44_DW.LastUAtTimeB = 0.0;
  helicopter44_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter44_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter44_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter44_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter44_M->Sizes.numY = (0);    /* Number of model outputs */
  helicopter44_M->Sizes.numU = (0);    /* Number of model inputs */
  helicopter44_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter44_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter44_M->Sizes.numBlocks = (61);/* Number of blocks */
  helicopter44_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helicopter44_M->Sizes.numBlockPrms = (155);/* Sum of parameter "widths" */
  return helicopter44_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
