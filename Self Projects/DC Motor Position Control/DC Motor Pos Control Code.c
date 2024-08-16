// Arduino Port Macros
#define ADC_CHANNEL_POT			(A0)
#define UNO_TO_IC_PIN1			(10U)
#define UNO_TO_IC_PIN2			(11U)

// Configurable motor-potentiometer parameters
#define MOTOR_MIN_ANGLE			(5U)
#define MOTOR_MAX_ANGLE			(350U)
#define MIN_POT_COUNTS			(0U)
#define MAX_POT_COUNTS			(1024U)
#define MOTOR_ROTATION_ANGLE	(180U)
#define MIN_PWM_COUNTS			(0U)
#define MAX_PWM_COUNTS			(255U)
#define NONLINEAR_START_ANGLE	(310U)
#define NONLINEAR_END_ANGLE		(330U)

// Configurable PID parameters
#define Kp						(1000.0)
#define Ki						(0.0)
#define Kd						(0.15)
#define DELAY_TIME_msec			(10U)

// Variables
float AdcCountsFromPot = 0.0;
float StartAngle = 0.0;
float EndAngle = 0.0;
float TargetAdcCounts = 0.0;
float CurrAdcCounts = 0.0;
float TimeInstance1_msec = 0.0, TimeInstance2_msec = 0.0, DeltaTime_sec = 0.0;
float ProportionalError = 0.0, prevProportionalError = 0.0;
float IntegralError = 0.0;
float DerivativeError;
float PidErrorVal = 0.0;
float test = 0;

// Functions
void Initialize(void);
void Stop_Motor(void);
void Rotate_Motor_Clockwise(float);
void Rotate_Motor_Anticlockwise(float);
float CountsToAngle(float);
float FindTargetAngle(float, float);
float AngleToCounts(float);

// Called only once
void setup()
{
	Initialize();
	AdcCountsFromPot = analogRead(ADC_CHANNEL_POT);
	StartAngle = CountsToAngle(AdcCountsFromPot);
	EndAngle = FindTargetAngle(StartAngle, MOTOR_ROTATION_ANGLE);
	TargetAdcCounts = AngleToCounts(EndAngle);
	prevProportionalError = TargetAdcCounts - AdcCountsFromPot;
	IntegralError = prevProportionalError;
	TimeInstance1_msec = millis();

	Serial.begin(9600);
}

// Called infinitely many times
void loop()
{
	CurrAdcCounts = analogRead(ADC_CHANNEL_POT);
	TimeInstance2_msec = millis();
	
	// Closed loop control logic (PID controller)
	DeltaTime_sec = (TimeInstance2_msec - TimeInstance1_msec)/1000.0;
	ProportionalError = TargetAdcCounts - CurrAdcCounts;
	IntegralError = IntegralError + ProportionalError*DeltaTime_sec;
	DerivativeError = (ProportionalError - prevProportionalError)/DeltaTime_sec;
	PidErrorVal = (Kp*ProportionalError) + (Ki*IntegralError) + (Kd*DerivativeError);
	test = PidErrorVal;

	if(PidErrorVal >= MAX_PWM_COUNTS)
	{
		PidErrorVal = MAX_PWM_COUNTS;
	}
	else if(PidErrorVal <= 0)
	{
		PidErrorVal = MIN_PWM_COUNTS;
	}
	else
	{
		/* No action here */
	}

	if(ProportionalError > 0)
	{
		Rotate_Motor_Clockwise(PidErrorVal);
	}
	else if(ProportionalError < 0)
	{
		Rotate_Motor_Anticlockwise(PidErrorVal);
	}
	else
	{
		Stop_Motor();
	}

	prevProportionalError = ProportionalError;
	TimeInstance1_msec = TimeInstance2_msec;
	delay(DELAY_TIME_msec);

	Serial.println("R");
	Serial.println(CurrAdcCounts);
	Serial.println(DeltaTime_sec);
	Serial.println(ProportionalError);
	Serial.println(IntegralError);
	Serial.println(DerivativeError);
	Serial.println(PidErrorVal);
	Serial.println(test);
}

void Initialize(void)
{
	pinMode(ADC_CHANNEL_POT, INPUT);
	pinMode(UNO_TO_IC_PIN1, OUTPUT);
	pinMode(UNO_TO_IC_PIN2, OUTPUT);
	Stop_Motor();
}

void Stop_Motor(void)
{
	analogWrite(UNO_TO_IC_PIN1, 0);
	analogWrite(UNO_TO_IC_PIN2, 0);
}

void Rotate_Motor_Clockwise(float dutycycle)
{
	analogWrite(UNO_TO_IC_PIN1, 0);
	analogWrite(UNO_TO_IC_PIN2, dutycycle);
}

void Rotate_Motor_Anticlockwise(float dutycycle)
{
	analogWrite(UNO_TO_IC_PIN1, dutycycle);
	analogWrite(UNO_TO_IC_PIN2, 0);
}

float CountsToAngle(float AdcCountsPot)
{
	return map(AdcCountsPot,MIN_POT_COUNTS,MAX_POT_COUNTS,MOTOR_MIN_ANGLE,MOTOR_MAX_ANGLE);
}

float FindTargetAngle(float StartTheta, float Delta)
{
	float TempVal;

	TempVal = StartTheta + Delta;
	if(TempVal >= 360)
	{
		TempVal = TempVal - 360;
	}
	else
	{
		/* No action here */ 
	}
	return TempVal;
}

float AngleToCounts(float EndTheta)
{
	return map(EndTheta,MOTOR_MIN_ANGLE,MOTOR_MAX_ANGLE,MIN_POT_COUNTS,MAX_POT_COUNTS);
}
