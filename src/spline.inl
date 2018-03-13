// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
    const T& position0, const T& position1, const T& tangent0,
    const T& tangent1, double normalizedTime, int derivative) {
	// TODO (Animation) Task 1a
	double t = normalizedTime;

	double h00, h10, h01, h11;

	if (derivative == 0) {
		h00 = 2. * t * t * t - 3. * t * t + 1.;
		h10 = t * t * t - 2. * t * t + t;
		h01 = -2. * t * t * t + 3. * t * t;
		h11 = t * t * t - t * t;
	}
	else if (derivative == 1) {
		h00 = 6. * t * t - 6. * t;
		h10 = 3. * t * t - 4. * t + 1;
		h01 = -6. * t * t + 6. * t;
		h11 = 3. * t * t - 2. * t;
	}
	else {
		h00 = 12. * t - 6.;
		h10 = 6. * t - 4.;
		h01 = -12. * t + 6.;
		h11 = 6. * t - 2.;
	}

	//cout << h00 << " / " << h10 << " / " << h01 << " / " << h11 << endl;

	return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

// Returns a state interpolated between the values directly before and after the
// given time.
template <class T>
inline T Spline<T>::evaluate(double time, int derivative) {
    // TODO (Animation) Task 1b
	//cout << "time: " << time << endl;
	// no knots in the spline
	if (knots.size() < 1) {
		//cout << "case 1" << endl;
		return T();
	}
	// only one knots in the spline
	else if (knots.size() == 1) {
		//cout << "case 2" << endl;
		if (derivative == 0) {
			return knots.begin()->second;
		}
		else {
			return T();
		}
	}
	// query time is less than intial knot
	else if (time <= knots.begin()->first) {
		//cout << "case 3" << endl;
		//cout << "time: " << time << endl;
		//cout << "first time: " << knots.begin()->first << endl;
		if (derivative == 0) {
			return knots.begin()->second;
		}
		else {
			return T();
		}
	}
	//query time is greater than the final knot
	else if (time >= prev(knots.end())->first) {
		//cout << "case 4" << endl;
		//cout << "time: " << time << endl;
		//cout << "last time: " << (--knots.end())->first << endl;
		if (derivative == 0) {
			return prev(knots.end())->second;
		}
		else {
			return T();
		}
	}
	//general case
	else {
		//cout << "case 5" << endl;
		KnotIter k0, k1, k2, k3;
		T p0, p1, p2, p3;
		T m1, m2;
		double norm_time;
		double t0, t1, t2, t3;

		k2 = knots.upper_bound(time);
		k1 = k2;
		k1--;
		t1 = k1->first;
		p1 = k1->second;
		t2 = k2->first;
		p2 = k2->second;

		norm_time = (time - t1) / (t2 - t1);

		// check if lower knot is the most left
		if (k1 == knots.begin()) {
			//cout << "case 5.1" << endl;
			p0 = p1 - (p2 - p1);
			t0 = t1 - (t2 - t1);
		}
		else {
			//cout << "case 5.2" << endl;
			k0 = k1;
			k0--;
			t0 = k0->first;
			p0 = k0->second;
		}
		
		// check if upper knot is the most right
		if (k2 == prev(knots.end())) {
			//cout << "case 5.3" << endl;
			p3 = p2 + (p2 - p1);
			t3 = t2 + (t2 - t1);
		}
		else {
			//cout << "case 5.4" << endl;
			k3 = k2;
			k3++;
			t3 = k3->first;
			p3 = k3->second;
		}

		double norm_t0 = (t0 - t1) / (t2 - t1); 
		double norm_t1 = 0;
		double norm_t2 = 1;
		double norm_t3 = (t3 - t1) / (t2 - t1); 


		m1 = (p2 - p0) / (norm_t2 - norm_t0);
		m2 = (p3 - p1) / (norm_t3 - norm_t1);

		return cubicSplineUnitInterval(p1, p2, m1, m2, norm_time, derivative);
	}
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance) {
  // Empty maps have no knots.
  if (knots.size() < 1) {
    return false;
  }

  // Look up the first element > or = to time.
  typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
  typename std::map<double, T>::iterator t1_iter;
  t1_iter = t2_iter;
  t1_iter--;

  if (t2_iter == knots.end()) {
    t2_iter = t1_iter;
  }

  // Handle tolerance bounds,
  // because we are working with floating point numbers.
  double t1 = (*t1_iter).first;
  double t2 = (*t2_iter).first;

  double d1 = fabs(t1 - time);
  double d2 = fabs(t2 - time);

  if (d1 < tolerance && d1 < d2) {
    knots.erase(t1_iter);
    return true;
  }

  if (d2 < tolerance && d2 < d1) {
    knots.erase(t2_iter);
    return t2;
  }

  return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue(double time, T value) {
  knots[time] = value;
}

template <class T>
inline T Spline<T>::operator()(double time) {
  return evaluate(time);
}
