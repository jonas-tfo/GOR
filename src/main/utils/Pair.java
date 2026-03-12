package main.utils;

public class Pair<T, V> {
	T v1;
	V v2;

	public Pair(T v1, V v2) {
		this.v1 = v1;
		this.v2 = v2;
	}

	public T getV1() {
		return this.v1;
	}

	public V getV2() {
		return this.v2;
	}

	public void setV1(T v) {
		this.v1 = v;
	}

	public void setV2(V v) {
		this.v2 = v;
	}

	@Override
	public String toString() {
		return "(" + this.v1.toString() + ", " + this.v2.toString() + ")";
	}
}
