# scurve library

S-curve motion profile generator library in C. Generates smooth motion trajectories with continuous acceleration and limited jerk. Features 7-segment profiles with independent acceleration/deceleration control. Suitable for embedded systems (STM32, etc) with no external dependencies. MIT licensed.

## SimpleFOC Example

```cpp
#include <Arduino.h>
#include <SimpleFOC.h>
#include "scurve.h"

// Configuration constants
#define PROFILE_UPDATE_INTERVAL_US 1000
#define PROFILE_UPDATE_INTERVAL_SEC (PROFILE_UPDATE_INTERVAL_US / 1000000.0f)

// Motor and driver setup (simplified)
BLDCMotor motor = BLDCMotor(7);
BLDCDriver3PWM driver = BLDCDriver3PWM(9, 5, 6, 8);
Encoder encoder = Encoder(2, 3, 8192);

// S-curve profile generator
unsigned char buffer[SCURVE_BUFFER_SIZE];
scurve_t* profile;

// Commander interface
Commander command = Commander(Serial);

// Command callbacks
void doTarget(char* cmd) {
    float target;
    command.scalar(&target, cmd);
    scurve_setPositionTarget(profile, target);
    if (scurve_startProfile(profile) != CURVE_SUCCESS) {
        Serial.println("Failed to start S-curve profile");
    }
}

void doAcceleration(char* cmd) {
    float acc;
    command.scalar(&acc, cmd);
    scurve_setAccelerationMax(profile, acc);
}

void doDeceleration(char* cmd) {
    float dec;
    command.scalar(&dec, cmd);
    scurve_setDecelerationMax(profile, dec);
}

void doJerk(char* cmd) {
    float j;
    command.scalar(&j, cmd);
    scurve_setJerkMax(profile, j);
}

void setup() {
    // Initialize motor and encoder (simplified)
    motor.linkSensor(&encoder);
    motor.linkDriver(&driver);
    motor.controller = MotionControlType::angle;
    motor.init();
    motor.initFOC();

    // Initialize S-curve profile generator
    profile = scurve_init(buffer, sizeof(buffer), PROFILE_UPDATE_INTERVAL_SEC);
    if (profile == NULL) {
        Serial.println("Failed to initialize S-curve profile");
        return;
    }

    // Configure S-curve parameters
    scurve_setPositionOutput(profile, 0.0f);
    scurve_setVelocityStart(profile, 0);
    scurve_setVelocityStop(profile, 0);
    scurve_setVelocityMax(profile, 250);       // rad/s
    scurve_setAccelerationMax(profile, 1600);  // rad/s²
    scurve_setDecelerationMax(profile, 1600);  // rad/s²
    scurve_setJerkMax(profile, 2500);          // rad/s³

    // Set up commander interface
    command.add('T', doTarget, "target angle");
    command.add('A', doAcceleration, "acceleration");
    command.add('D', doDeceleration, "deceleration");
    command.add('J', doJerk, "jerk");
    
    Serial.begin(115200);
}

uint32_t last_profile_update_us = 0;

void loop() {
    // Update S-curve profile at fixed interval
    if (micros() - last_profile_update_us > PROFILE_UPDATE_INTERVAL_US) {
        scurve_state_e state = scurve_run(profile);
        if (state == CURVE_BUSY || state == CURVE_ONEND) {
            motor.target = scurve_getPositionOutput(profile);
        }
        last_profile_update_us = micros();
    }

    // Run motor control
    motor.loopFOC();
    motor.move();

    // Handle user commands
    command.run();
}
```

## License

MIT License

Copyright (c) 2025 