<template>
  <div class="app">
    <header>
      <img alt="Vue logo" class="logo" src="./assets/logo.svg" width="125" height="125" />
      <div class="wrapper">
        <h1>Genomic Risk Analysis</h1>
      </div>
    </header>

    <main>
      <RiskViewer />
    </main>
  </div>
</template>

<script setup lang="ts">
import RiskViewer from './components/RiskViewer.vue'
</script>

<style>
:root {
  --primary-color: #1976d2;
  --secondary-color: #2c3e50;
  --background-color: #f8f9fa;
  --text-color: #2c3e50;
  --border-radius: 8px;
  --box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

body {
  margin: 0;
  padding: 0;
  background-color: var(--background-color);
  color: var(--text-color);
  font-family: 'Arial', sans-serif;
}

.app {
  min-height: 100vh;
  display: flex;
  flex-direction: column;
}

header {
  background-color: white;
  padding: 2rem;
  box-shadow: var(--box-shadow);
  margin-bottom: 2rem;
}

.logo {
  display: block;
  margin: 0 auto 2rem;
  transition: transform 0.3s ease;
}

.logo:hover {
  transform: scale(1.05);
}

.wrapper {
  text-align: center;
}

h1 {
  color: var(--secondary-color);
  margin: 0;
  font-size: 2rem;
}

main {
  flex: 1;
  padding: 0 1rem;
  max-width: 1200px;
  margin: 0 auto;
  width: 100%;
}

@media (min-width: 1024px) {
  header {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 2rem;
  }

  .logo {
    margin: 0;
  }

  .wrapper {
    text-align: left;
  }
}
</style> 