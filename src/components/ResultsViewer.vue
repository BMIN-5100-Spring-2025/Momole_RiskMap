<template>
  <div class="results-viewer">
    <h2>Genomic Risk Analysis Results</h2>
    <div v-if="loading" class="loading">
      <div class="spinner"></div>
      <p>Loading results...</p>
    </div>
    <div v-else-if="error" class="error">
      <div class="error-icon">⚠️</div>
      <p>{{ error }}</p>
    </div>
    <div v-else class="results">
      <div class="patient-info">
        <div class="info-item">
          <span class="label">Patient ID:</span>
          <span class="value">{{ parsedResults.patient_id }}</span>
        </div>
        <div class="info-item">
          <span class="label">Analysis Date:</span>
          <span class="value">{{ parsedResults.analysis_date }}</span>
        </div>
      </div>

      <div class="risk-factors">
        <h3>Risk Factors</h3>
        <div class="risk-factor-list">
          <div v-for="(factor, index) in parsedResults.risk_factors" :key="index" 
               :class="['risk-factor', factor.risk_level.toLowerCase()]">
            <div class="risk-header">
              <span class="gene">{{ factor.gene }}</span>
              <span class="variant">{{ factor.variant }}</span>
            </div>
            <div class="risk-details">
              <span class="risk-level">{{ factor.risk_level }} Risk</span>
              <span class="confidence">Confidence: {{ (factor.confidence * 100).toFixed(1) }}%</span>
            </div>
          </div>
        </div>
      </div>

      <div class="summary">
        <h3>Summary</h3>
        <div class="summary-stats">
          <div class="stat-item">
            <span class="stat-value">{{ parsedResults.summary.total_risk_factors }}</span>
            <span class="stat-label">Total Risk Factors</span>
          </div>
          <div class="stat-item high">
            <span class="stat-value">{{ parsedResults.summary.high_risk_count }}</span>
            <span class="stat-label">High Risk</span>
          </div>
          <div class="stat-item medium">
            <span class="stat-value">{{ parsedResults.summary.medium_risk_count }}</span>
            <span class="stat-label">Medium Risk</span>
          </div>
          <div class="stat-item low">
            <span class="stat-value">{{ parsedResults.summary.low_risk_count }}</span>
            <span class="stat-label">Low Risk</span>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, computed } from 'vue'

const results = ref('')
const loading = ref(true)
const error = ref('')

const parsedResults = computed(() => {
  try {
    return JSON.parse(results.value)
  } catch {
    return {}
  }
})

onMounted(async () => {
  try {
    const response = await fetch('/results.json')
    if (!response.ok) throw new Error('Failed to load results')
    results.value = await response.text()
  } catch (err) {
    error.value = err instanceof Error ? err.message : 'An error occurred'
  } finally {
    loading.value = false
  }
})
</script>

<style scoped>
.results-viewer {
  padding: 2rem;
  max-width: 1000px;
  margin: 0 auto;
  font-family: 'Arial', sans-serif;
}

h2 {
  color: #2c3e50;
  text-align: center;
  margin-bottom: 2rem;
  font-size: 2rem;
}

h3 {
  color: #2c3e50;
  margin: 1.5rem 0 1rem;
  font-size: 1.5rem;
}

.loading {
  text-align: center;
  padding: 2rem;
}

.spinner {
  width: 40px;
  height: 40px;
  border: 4px solid #f3f3f3;
  border-top: 4px solid #3498db;
  border-radius: 50%;
  margin: 0 auto 1rem;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

.error {
  text-align: center;
  padding: 2rem;
  background-color: #ffebee;
  border-radius: 8px;
  color: #c62828;
}

.error-icon {
  font-size: 2rem;
  margin-bottom: 1rem;
}

.patient-info {
  background-color: #e3f2fd;
  padding: 1.5rem;
  border-radius: 8px;
  margin-bottom: 2rem;
  display: flex;
  gap: 2rem;
}

.info-item {
  display: flex;
  flex-direction: column;
}

.label {
  font-weight: bold;
  color: #1976d2;
  margin-bottom: 0.5rem;
}

.value {
  color: #2c3e50;
}

.risk-factors {
  margin-bottom: 2rem;
}

.risk-factor-list {
  display: grid;
  gap: 1rem;
}

.risk-factor {
  background-color: white;
  border-radius: 8px;
  padding: 1.5rem;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.risk-factor.high {
  border-left: 4px solid #e53935;
}

.risk-factor.medium {
  border-left: 4px solid #fb8c00;
}

.risk-factor.low {
  border-left: 4px solid #43a047;
}

.risk-header {
  display: flex;
  justify-content: space-between;
  margin-bottom: 1rem;
}

.gene {
  font-weight: bold;
  font-size: 1.2rem;
  color: #2c3e50;
}

.variant {
  color: #666;
  font-family: monospace;
}

.risk-details {
  display: flex;
  justify-content: space-between;
  color: #666;
}

.risk-level {
  font-weight: bold;
}

.summary {
  background-color: #f5f5f5;
  padding: 1.5rem;
  border-radius: 8px;
}

.summary-stats {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 1rem;
  margin-top: 1rem;
}

.stat-item {
  background-color: white;
  padding: 1rem;
  border-radius: 8px;
  text-align: center;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.stat-item.high {
  border-top: 4px solid #e53935;
}

.stat-item.medium {
  border-top: 4px solid #fb8c00;
}

.stat-item.low {
  border-top: 4px solid #43a047;
}

.stat-value {
  display: block;
  font-size: 2rem;
  font-weight: bold;
  color: #2c3e50;
  margin-bottom: 0.5rem;
}

.stat-label {
  color: #666;
  font-size: 0.9rem;
}
</style> 