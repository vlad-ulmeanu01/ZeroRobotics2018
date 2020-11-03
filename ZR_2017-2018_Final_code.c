// Macros (compiler inserts the value where we write the key)

// - States
#define GO_TO_SAMPLE 0
#define COLLECT_SAMPLE 1
#define RETURN_SAMPLES 2

// - Game data
#define WIDTH 16
#define HEIGHT 20

// - Drilling
#define MAXIMUM_NUMBER_OF_ITEMS_WHICH_CAN_BE_HELD 5
#define NUMBER_OF_ITEMS_TO_RETURN_AT 4
#define NUMBER_OF_DRILLS_TO_STOP_DRILLING_A_SQUARE_AT 3

// - Movement
#define SET_POSITION_TARGET_THRESHOLD 0.1f
#define MOVE_SPEED_COEFFICIENT 0.21f
#define SPHERE_DROP_RADIUS 0.2f

// - Other
#define SQUARE_SIDE_LENGTH 0.08f
#define FIRST_SHELL_RADIUS 0.13f
#define SECOND_SHELL_RADIUS 0.23f
#define FIRST_SHELL_DISTANCE_DIVISOR 3
#define SECOND_SHELL_DISTANCE_DIVISOR 3

float myZRState[12];
float otherZRState[12];
float positionToDrillAt[3];
float positionToDrillAtEnemy[3];
float vector[3]; // A multi-purpose vector.
float positionToDropSamplesAt[3];
float up[3];
float onXYPlaneAttitude[3];
float bestSquareFoundSoFar[2];
float coordinatesOfHeldSamples[MAXIMUM_NUMBER_OF_ITEMS_WHICH_CAN_BE_HELD][2];
float bestConcentrationFoundSoFar;

int step;

void init() {
    up[0] =   0.0;
    up[1] =   0.0;
    up[2] = - 1.0;

    step = 0;

    // So not squares will be close enough to begin with.
    bestSquareFoundSoFar[0] = 10.f;
    bestSquareFoundSoFar[1] = 10.f;

    // Better than 0.1, but worse than 0.3.
    bestConcentrationFoundSoFar = 0.2f;

    DEBUG ( ( "cheeki breeki" ) );
}

void loop() {
    // DEBUG(
    //     (
    //         "bestSquareFoundSoFar: [%f, %f]",
    //         bestSquareFoundSoFar[0],
    //         bestSquareFoundSoFar[1])
    // );
    // DEBUG(
    //     ("bestConcentrationFoundSoFar: %f", bestConcentrationFoundSoFar)
    // );

    api.getMyZRState(myZRState);
    onXYPlaneAttitude[0] = myZRState[6];
    onXYPlaneAttitude[1] = myZRState[7];
    onXYPlaneAttitude[2] = 0.0f;
    mathVecNormalize(onXYPlaneAttitude, 3);

    //DEBUG(("step: %d", step));

    if (game.getFuelRemaining () > 0.02f || game.getNumSamplesHeld() > 0) {
        if (step == GO_TO_SAMPLE) {
            setPositionToDrillAt();
            move(positionToDrillAt);
            api.setAttitudeTarget(onXYPlaneAttitude);

            if (canDrill()) {
                DEBUG ( ( "orientation %f", myZRState[8] ) );
                step = COLLECT_SAMPLE;
            }

        }

        if (step == COLLECT_SAMPLE) {
            move(positionToDrillAt);

            if (!game.getDrillEnabled()) {
                game.startDrill();
            }

            float desiredAttitude[3];
            mathVecCross(desiredAttitude, up, onXYPlaneAttitude);
            mathVecNormalize(desiredAttitude, 3);
            api.setAttitudeTarget(desiredAttitude);

            if (game.checkSample()) {
                if (
                        game.getNumSamplesHeld()
                            >= MAXIMUM_NUMBER_OF_ITEMS_WHICH_CAN_BE_HELD
                ) {
                    game.dropSample(0);
                }

                int sampleIndex = game.pickupSample();
                coordinatesOfHeldSamples[sampleIndex][0] = positionToDrillAt[0];
                coordinatesOfHeldSamples[sampleIndex][1] = positionToDrillAt[1];
            }


            if (
                    game.isGeyserHere (myZRState)
                    || game.getDrills (myZRState)
                        >= NUMBER_OF_DRILLS_TO_STOP_DRILLING_A_SQUARE_AT
            ) {
                game.stopDrill();
                setPositionToDrillAt();

                if (game.getNumSamplesHeld() >= NUMBER_OF_ITEMS_TO_RETURN_AT) {
                    step = RETURN_SAMPLES;
                }

                else {
                    step = GO_TO_SAMPLE;
                }
            }

        }

        if (step == RETURN_SAMPLES) {
            float distanceFromOrigin = mathVecMagnitude(myZRState, 3);

            for (int index = 0; index < 3; index ++) {
                positionToDropSamplesAt[index] = myZRState[index]
                    * SPHERE_DROP_RADIUS
                    / distanceFromOrigin;
            }

            api.setAttitudeTarget(up);

            move(positionToDropSamplesAt);

            if (game.atBaseStation()) {
                for (int sampleIndex = 0; sampleIndex < 5; sampleIndex++) {
                    float concentration = game.dropSample(sampleIndex);

                    if (concentration > bestConcentrationFoundSoFar) {
                        bestSquareFoundSoFar[0] =
                            coordinatesOfHeldSamples[sampleIndex][0];
                        bestSquareFoundSoFar[1] =
                            coordinatesOfHeldSamples[sampleIndex][1];
                        bestConcentrationFoundSoFar = concentration;
                    }
                }
                step = GO_TO_SAMPLE;
            }
        }
    }

    else  {
        // Going back to the base, as otherwise we will run out of fuel and
        // crash.
        game.stopDrill ();
        float origin[3] = {0.0, 0.0, 0.0};
        move(origin);
    }

    if (game.getDrillError()) {
        game.stopDrill();
        step = GO_TO_SAMPLE;
    }

    if (shouldReturnToBaseBecauseOfTime() == true && step != RETURN_SAMPLES) {
        game.stopDrill();
        step = RETURN_SAMPLES;
    }
}

bool isGeyserNear(float positionToCheckAround[3]) {
    // Returns true if the square or any square arround it has an active
    // geyser.
    float positionBeingChecked[2];

    for (int xOffset = -1; xOffset <= 1; xOffset++) {
        for (int yOffset = -1; yOffset <= 1; yOffset++) {
            positionBeingChecked[0] = positionToCheckAround[0] +
                xOffset * SQUARE_SIDE_LENGTH;
            positionBeingChecked[1] = positionToCheckAround[1] +
                yOffset * SQUARE_SIDE_LENGTH;

            if (game.isGeyserHere(positionBeingChecked)) {
                return true;
            }
        }
    }

    // None of them did.
    return false;
}

void setPositionToDrillAt ( ) {
    float smallestDistanceSoFar = 10.0f;
    // There should certainly be something closer than that.
    float coordinatesBeingTested[3];
    // Sphere radius = 0.11, and the distance we want to be above the surface =
    // 0.02 (to be right in the middle of the 0.04 range), and we are only
    // looking at samples of a height of 0.4, so the z component is always
    // 0.4 - 0.11 - 0.2 = 0.27.
    coordinatesBeingTested[2] = 0.2675f; /// this should no longer allow hitting the ground

    for ( int i = 0; i < WIDTH; i++ )
    {
        for ( int j = 0; j < HEIGHT; j++ )
        {
            coordinatesBeingTested[0] = -0.6 + 0.08 * i;
            coordinatesBeingTested[1] = -0.76 + 0.08 * j;

            if (6 <= i && i <= 9 && 8 <= j && j <= 11) {
                // These coordinates are under the base station
                continue;
            }

            api.getMyZRState ( myZRState );
            api.getOtherZRState ( otherZRState );

            mathVecSubtract ( vector, coordinatesBeingTested, otherZRState, 3 );
            float distToTheOtherSatelliteWhenDrilling = mathVecMagnitude ( vector, 3 );

            if ( distToTheOtherSatelliteWhenDrilling < 0.31f && coordinatesBeingTested[2] >= otherZRState[2] - 0.015f  )
            {
                //DEBUG ( ( "oth sq %f delta %f", distToTheOtherSatelliteWhenDrilling, coordinatesBeingTested[2] - otherZRState[2] + 0.015f ) );
                continue;
            }

            // Set vector to be the displacement of coordinatesBeingTested
            // from myZRState.
            mathVecSubtract(vector, coordinatesBeingTested, myZRState, 3);
            // The distance is sum of the distance from us to the square, and
            // the distance from the origin to the square.
            float distanceForThisSquare = mathVecMagnitude(vector, 3) +
                mathVecMagnitude(coordinatesBeingTested, 3);

            float distanceToBestSquareFoundSoFar =
                distanceOfPointToBestSquareFoundSoFar(coordinatesBeingTested);

            if (distanceToBestSquareFoundSoFar <= SECOND_SHELL_RADIUS) {
                // If this square is in the second shell of the best square we
                // have found so far, pretend it is much closer.
                distanceForThisSquare /= FIRST_SHELL_DISTANCE_DIVISOR;
            }
            if (distanceToBestSquareFoundSoFar <= FIRST_SHELL_RADIUS) {
                // If it is in the first shell, pretend it is even closer again
                // (if this is happening, the previous if statement also
                // happened).
                distanceForThisSquare /= SECOND_SHELL_DISTANCE_DIVISOR;
            }

            if (
                game.getTerrainHeight(coordinatesBeingTested) == 0.4f
                && distanceForThisSquare < smallestDistanceSoFar
                && !isGeyserNear(coordinatesBeingTested)
                && game.getDrills(coordinatesBeingTested) == 0
            ) {
                smallestDistanceSoFar = distanceForThisSquare;
                for ( int index = 0; index < 3; index++ )
                {
                    positionToDrillAt[index] = coordinatesBeingTested[index];
                }
            }
        }
    }
}

float distanceOfPointToBestSquareFoundSoFar(float point[2]) {
    mathVecSubtract(
        vector,
        point,
        bestSquareFoundSoFar,
        2
    );
    float distanceToBestSquareFoundSoFar = mathVecMagnitude(
        vector,
        2
    );

    return distanceToBestSquareFoundSoFar;
}

bool canDrill() {
    // If the z component of the sphere's attitude is smaller than 0.19, the
    // the angle that attitude makes with the xy plane is smaller than 11
    // degrees. You can see this by drawing a right angled triangle with the
    // attitude as the hypotenuse, the xy plane as the base, and hence the z
    // component of attitude as the remaining side.

    api.getMyZRState(myZRState);

    DEBUG ( ( "%d %d %d %d %d %d",
    myZRState[2] <= 0.288f && myZRState[2] >= 0.252f,
    fabsf(myZRState[0] - positionToDrillAt[0]) <= 0.04f,
    fabsf(myZRState[1] - positionToDrillAt[1]) <= 0.04f,
    mathVecMagnitude (&myZRState[9], 3) < 0.04f,
    fabsf (myZRState[8]) < 0.18f,
    mathVecMagnitude (&myZRState[3], 3) < 0.01f ) );

    return (
        // The sattelite lines within the allowed drilling rectangle.
        //fabsf(myZRState[2] - positionToDrillAt[2]) <= 0.019f
        myZRState[2] <= 0.288f && myZRState[2] >= 0.252f
        && fabsf(myZRState[0] - positionToDrillAt[0]) <= 0.04f
        && fabsf(myZRState[1] - positionToDrillAt[1]) <= 0.04f
        // The rotational speed is within restrictions.
        && mathVecMagnitude (&myZRState[9], 3) < 0.04f
        // The angle with the xy plane is within restrictions.
        && fabsf (myZRState[8]) < 0.18f
        // The speed is within restrictions.
        && mathVecMagnitude (&myZRState[3], 3) < 0.01f
    );
}

void normal_move ( float positionToMoveTo[3], float move_coeff )
{
    float between[3];

    api.getMyZRState(myZRState);

    for ( int index = 0; index < 3; index++ )
    {
        between[index] = move_coeff * (
            positionToMoveTo[index] - myZRState[index]
        );
    }

    api.setVelocityTarget (between);
}

void move ( float positionToMoveTo[3] ) {
    // Set vector to be the displacement of positionToMove to from myZRState.
    mathVecSubtract(vector, positionToMoveTo, myZRState, 3);

    api.getMyZRState(myZRState);

    if (game.isGeyserHere(myZRState))
    {
        DEBUG(("HERE"));
        // Flatten the vector to the xy plane, such that you don't waste fuel
        // pushing yourself in the direction the geyser is already pushing you.
        vector[2] = 0.f;
        // Normalise the force such that it is larger than you could possible
        // apply.
        mathVecNormalize(vector, 3);
        api.setForces(vector);
    }
    else if (mathVecMagnitude(vector, 3) > SET_POSITION_TARGET_THRESHOLD )
    {
        normal_move ( positionToMoveTo, MOVE_SPEED_COEFFICIENT );
    }
    else
    {
        api.setPositionTarget(positionToMoveTo);
    }
}

bool shouldReturnToBaseBecauseOfTime() {
    api.getMyZRState(myZRState);

    float distanceOnXYPlane = mathVecMagnitude(myZRState, 2);

    if (
            ( game.getNumSamplesHeld() >= 3 )
            && (
                game.getFuelRemaining () < 0.15
                || (
                    api.getTime () >= 160
                    && distanceOnXYPlane < 8 * 0.08
                )
                || (
                    api.getTime () >= 155
                    && distanceOnXYPlane >= 8 * 0.08
                )
            )
    )
    {
        DEBUG (("Returning to base because of time."));
        return true;
    }
    else
        return false;
}
