<template>
  <div class="risk-viewer">
    <div v-if="loading" class="loading">
      <div class="spinner"></div>
      <p>Loading risk analysis...</p>
    </div>
    
    <div v-else-if="error" class="error">
      <p class="error-message">{{ error }}</p>
      <button @click="retry" class="retry-button">Retry</button>
    </div>
    
    <div v-else-if="!hasData" class="no-data">
      <p>No risk analysis data available.</p>
      <button @click="retry" class="retry-button">Refresh</button>
    </div>
    
    <div v-else class="results">
      <div class="patient-info">
        <h2>Patient Information</h2>
        <div class="info-grid">
          <div class="info-item">
            <span class="label">Patient ID:</span>
            <span class="value">{{ patientId }}</span>
          </div>
          <div class="info-item">
            <span class="label">Analysis Date:</span>
            <span class="value">{{ analysisDate }}</span>
          </div>
        </div>
      </div>
      
      <div class="risk-factors">
        <h2>Risk Factors</h2>
        <div class="risk-grid">
          <div v-for="(factor, index) in riskFactors" :key="index" 
               class="risk-item" :class="factor.risk_level.toLowerCase()">
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
        <h2>Summary</h2>
        <div class="summary-stats">
          <div class="stat-item">
            <span class="stat-value">{{ summary.total_risk_factors }}</span>
            <span class="stat-label">Total Risk Factors</span>
          </div>
          <div class="stat-item high">
            <span class="stat-value">{{ summary.high_risk_count }}</span>
            <span class="stat-label">High Risk</span>
          </div>
          <div class="stat-item medium">
            <span class="stat-value">{{ summary.medium_risk_count }}</span>
            <span class="stat-label">Medium Risk</span>
          </div>
          <div class="stat-item low">
            <span class="stat-value">{{ summary.low_risk_count }}</span>
            <span class="stat-label">Low Risk</span>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, computed } from 'vue'

const loading = ref(true)
const error = ref(null)
const patientId = ref('')
const analysisDate = ref('')
const riskFactors = ref([])
const summary = ref({
  total_risk_factors: 0,
  high_risk_count: 0,
  medium_risk_count: 0,
  low_risk_count: 0
})

const hasData = computed(() => {
  return riskFactors.value.length > 0
})

const fetchData = async () => {
  try {
    loading.value = true
    error.value = null
    const response = await fetch('http://localhost:5000/api/risk-analysis')
    
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`)
    }
    
    const data = await response.json()
    
    if (data.error) {
      throw new Error(data.error)
    }
    
    patientId.value = data.patient_id
    analysisDate.value = data.analysis_date
    riskFactors.value = data.risk_factors
    summary.value = data.summary
  } catch (err) {
    error.value = 'Failed to load risk analysis data. Please try again.'
    console.error('Error fetching data:', err)
  } finally {
    loading.value = false
  }
}

const retry = () => {
  fetchData()
}

// Initial fetch
fetchData()

// Set up auto-refresh every 5 minutes
setInterval(fetchData, 5 * 60 * 1000)
</script>

<style scoped>
.risk-viewer {
  max-width: 1200px;
  margin: 0 auto;
  padding: 2rem;
}

.loading {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  min-height: 400px;
}

.spinner {
  width: 50px;
  height: 50px;
  border: 5px solid #f3f3f3;
  border-top: 5px solid #3498db;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

.error, .no-data {
  text-align: center;
  padding: 2rem;
}

.error-message {
  color: #e74c3c;
  margin-bottom: 1rem;
}

.retry-button {
  padding: 0.5rem 1rem;
  background-color: #3498db;
  color: white;
  border: none;
  border-radius: 4px;
  cursor: pointer;
  transition: background-color 0.3s ease;
}

.retry-button:hover {
  background-color: #2980b9;
}

.patient-info, .risk-factors, .summary {
  margin-bottom: 2rem;
  background-color: white;
  padding: 1.5rem;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

h2 {
  color: #2c3e50;
  margin-bottom: 1rem;
  font-size: 1.5rem;
}

.info-grid, .risk-grid {
  display: grid;
  gap: 1rem;
}

.info-item {
  display: flex;
  justify-content: space-between;
  padding: 0.5rem;
  background-color: #f8f9fa;
  border-radius: 4px;
}

.label {
  font-weight: bold;
  color: #2c3e50;
}

.value {
  color: #34495e;
}

.risk-item {
  padding: 1rem;
  background-color: white;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  transition: transform 0.2s ease;
}

.risk-item:hover {
  transform: translateY(-2px);
}

.risk-item.high {
  border-left: 4px solid #e53935;
}

.risk-item.medium {
  border-left: 4px solid #fb8c00;
}

.risk-item.low {
  border-left: 4px solid #43a047;
}

.risk-header {
  display: flex;
  justify-content: space-between;
  margin-bottom: 0.5rem;
}

.gene {
  font-weight: bold;
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

.summary-stats {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 1rem;
}

.stat-item {
  background-color: white;
  padding: 1rem;
  border-radius: 8px;
  text-align: center;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  transition: transform 0.2s ease;
}

.stat-item:hover {
  transform: translateY(-2px);
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

@media (max-width: 768px) {
  .risk-viewer {
    padding: 1rem;
  }
  
  .summary-stats {
    grid-template-columns: 1fr;
  }
  
  .risk-item, .stat-item {
    padding: 0.75rem;
  }
}
</style> 